[CmdletBinding()]
param(
    [string]$SourceDirectory,
    [string]$RunDirectory,
    [ValidateSet('L','A')][string]$Stage='L',
    [switch]$PrepareOnly
)

$ErrorActionPreference='Stop'
if (-not $SourceDirectory) {
    $SourceDirectory = Join-Path (Split-Path -Parent $PSScriptRoot) 'fortran'
}
$env:PSModulePath=(Join-Path $PSHOME 'Modules')+[IO.Path]::PathSeparator+$env:PSModulePath
Import-Module Microsoft.PowerShell.Management
Import-Module Microsoft.PowerShell.Utility
$utf8=New-Object Text.UTF8Encoding($false)
$source=[IO.Path]::GetFullPath($SourceDirectory)
function Read-ControlRows([string]$path,[int]$count) {
    $rows=@([IO.File]::ReadAllLines($path)|ForEach-Object{($_ -split '[!#]',2)[0].Trim()}|Where-Object{$_})
    if($rows.Count -lt $count){throw "Missing control rows: $path"}
    return ,$rows[0..($count-1)]
}
function Hash-File([string]$path) {
    $stream=[IO.File]::Open($path,'Open','Read',([IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete))
    $sha=[Security.Cryptography.SHA256]::Create()
    try{return ([BitConverter]::ToString($sha.ComputeHash($stream))).Replace('-','')}finally{$stream.Dispose();$sha.Dispose()}
}
function Save-State([string]$path,$value) {
    $temp=$path+'.tmp-'+[guid]::NewGuid().ToString('N')
    [IO.File]::WriteAllText($temp,($value|ConvertTo-Json -Depth 12),$utf8)
    if([IO.File]::Exists($path)){
        $backup=$path+'.bak-'+[guid]::NewGuid().ToString('N')
        try{[IO.File]::Replace($temp,$path,$backup)}finally{if([IO.File]::Exists($backup)){[IO.File]::Delete($backup)}}
    }else{[IO.File]::Move($temp,$path)}
}
function Read-State([string]$path) {
    $stream=[IO.File]::Open($path,'Open','Read',([IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete))
    $reader=New-Object IO.StreamReader($stream,[Text.Encoding]::UTF8,$true)
    try{return ($reader.ReadToEnd()|ConvertFrom-Json)}finally{$reader.Dispose();$stream.Dispose()}
}

$control=Read-ControlRows (Join-Path $source 'Input file/INPUT.in') 7
$nfile=[int]$control[0];$eventsPerTask=[int]$control[1];$name=$control[2]
if($nfile -lt 1 -or $nfile -gt 5 -or $eventsPerTask -lt 1){throw 'This five-table batch requires NFILE=1..5 and positive NRUN_PER, giving at most 25 calculations.'}
if($name -in @('.','..') -or $name.IndexOfAny([IO.Path]::GetInvalidFileNameChars()) -ge 0){throw 'Invalid NAME in INPUT.in.'}
if([int]$control[6] -ne 0){throw 'FAIL must be zero; incomplete physics events are not counted.'}
if(-not $RunDirectory){$RunDirectory=Join-Path $source ('results/'+$name+'_mass_tables')}
$run=[IO.Path]::GetFullPath($RunDirectory)
$binary=Join-Path $source 'x64/Release/win_batch'
$files=New-Object 'System.Collections.Generic.List[IO.FileInfo]'
foreach($folder in @('Input file','gdat','tbl')){
    Get-ChildItem -LiteralPath (Join-Path $source $folder) -Recurse -File |ForEach-Object{$files.Add($_)}
}
foreach($file in @('pini.dat','tini.dat')){$files.Add((Get-Item -LiteralPath (Join-Path $source $file)))}
foreach($file in @('LQMD.exe','exfrag.exe','ctrlL_win.exe','ctrlA_win.exe')){$files.Add((Get-Item -LiteralPath (Join-Path $binary $file)))}
$fingerprint=[ordered]@{}
foreach($file in $files|Sort-Object FullName){$fingerprint[$file.FullName.Substring($source.TrimEnd('\').Length+1)]=Hash-File $file.FullName}
foreach($file in @('windows_event_controller.ps1','windows_controller_results.ps1')){$fingerprint['scripts/'+$file]=Hash-File (Join-Path $PSScriptRoot $file)}
$identity=[ordered]@{source=$source;events_per_model=$nfile*$eventsPerTask;tasks_per_model=$nfile;models=@(1,2,3,4,5);input_hashes=$fingerprint}
$pathBytes=[Text.Encoding]::UTF8.GetBytes($run.ToLowerInvariant())
$sha=[Security.Cryptography.SHA256]::Create()
try{$suffix=([BitConverter]::ToString($sha.ComputeHash($pathBytes))).Replace('-','').Substring(0,24)}finally{$sha.Dispose()}
$mutex=New-Object Threading.Mutex($false,('Local\SpallationFiveTableBatch_'+$suffix))
$locked=$false
try{
    try{$locked=$mutex.WaitOne(0)}catch [Threading.AbandonedMutexException]{$locked=$true}
    if(-not $locked){throw 'This five-table batch is already controlled by another process.'}
    $manifestPath=Join-Path $run 'batch-manifest.json'
    if(Test-Path -LiteralPath $manifestPath){
        $old=Read-State $manifestPath
        if(($old.identity|ConvertTo-Json -Depth 12 -Compress) -cne ($identity|ConvertTo-Json -Depth 12 -Compress)){throw 'Source controls, templates or executables changed; use a new NAME or RunDirectory.'}
    }else{
        if((Test-Path -LiteralPath $run) -and @(Get-ChildItem -LiteralPath $run -Force).Count){throw 'RunDirectory is not empty and is not this batch.'}
        New-Item -ItemType Directory -Path $run -Force|Out-Null
        Save-State $manifestPath ([ordered]@{identity=$identity;created_utc=[DateTime]::UtcNow.ToString('o');notice='Five independent LQMD ensembles; framework 1 still freezes dynamics to M0.'})
    }
    foreach($model in 1..5){
        $modelRoot=Join-Path $run ('model-'+$model)
        $stamp=Join-Path $modelRoot 'prepared.json'
        if(-not(Test-Path -LiteralPath $stamp)){
            if((Test-Path -LiteralPath $modelRoot) -and @(Get-ChildItem -LiteralPath $modelRoot -Force).Count){throw "Incomplete existing model workspace: $modelRoot"}
            New-Item -ItemType Directory -Path $modelRoot -Force|Out-Null
            foreach($folder in @('Input file','gdat','tbl')){Copy-Item -LiteralPath (Join-Path $source $folder) -Destination $modelRoot -Recurse}
            foreach($file in @('pini.dat','tini.dat')){Copy-Item -LiteralPath (Join-Path $source $file) -Destination $modelRoot}
            foreach($file in @('LQMD.exe','exfrag.exe','ctrlL_win.exe','ctrlA_win.exe')){Copy-Item -LiteralPath (Join-Path $binary $file) -Destination $modelRoot}
            $scriptRoot=Join-Path $modelRoot 'scripts'
            New-Item -ItemType Directory -Path $scriptRoot|Out-Null
            foreach($file in @('windows_event_controller.ps1','windows_controller_results.ps1')){Copy-Item -LiteralPath (Join-Path $PSScriptRoot $file) -Destination $scriptRoot}
            $selector=Join-Path $modelRoot 'Input file/m_select.in'
            $lines=[IO.File]::ReadAllLines($selector);$row=0
            for($i=0;$i -lt $lines.Count;$i++){
                if($lines[$i] -match '^\s*\d+\s*(?:!.*)?$'){
                    $row++;if($row -eq 1){$lines[$i]=[regex]::Replace($lines[$i],'^\s*\d+',[string]$model)}
                    if($row -eq 2){$lines[$i]=[regex]::Replace($lines[$i],'^\s*\d+','1')}
                    if($row -eq 3){break}
                }
            }
            if($row -lt 3){throw 'm_select.in requires three numeric values.'}
            [IO.File]::WriteAllLines($selector,$lines,$utf8)
            Save-State $stamp ([ordered]@{model=$model;events=$nfile*$eventsPerTask;tasks=$nfile;mass_framework=1;created_utc=[DateTime]::UtcNow.ToString('o')})
        }
    }
    $statePath=Join-Path $run 'batch-status.json'
    $state=[ordered]@{stage='prepared';pid=$PID;controller_stage=$Stage;events_per_model=$nfile*$eventsPerTask;tasks_per_model=$nfile;total_compute_limit=5*$nfile;models=@();updated_utc=[DateTime]::UtcNow.ToString('o')}
    Save-State $statePath $state
    if($PrepareOnly){Write-Output "Prepared five table workspaces at $run";return}
    $controllers=New-Object 'System.Collections.Generic.List[object]'
    foreach($model in 1..5){
        $modelRoot=Join-Path $run ('model-'+$model)
        $out=Join-Path $modelRoot ('controller-'+[DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfff')+'.stdout.log')
        $err=$out.Replace('.stdout.log','.stderr.log')
        $controller=if($Stage -eq 'L'){'ctrlL_win.exe'}else{'ctrlA_win.exe'}
        $process=Start-Process -FilePath (Join-Path $modelRoot $controller) -WorkingDirectory $modelRoot -PassThru -WindowStyle Hidden -RedirectStandardOutput $out -RedirectStandardError $err
        $null=$process.Handle
        $controllers.Add([pscustomobject]@{Model=$model;Process=$process;Output=$out;Error=$err;Root=$modelRoot})
        $state.models+=@{model=$model;controller_pid=$process.Id;stdout=$out;stderr=$err;status='running'}
    }
    $state.stage='running';Save-State $statePath $state
    Write-Output "Started five stage $Stage controllers, each NFILE=$nfile, NRUN_PER=$eventsPerTask; at most $($nfile*5) LQMD/exfrag calculations."
    while($true){
        $active=0;$failed=@()
        foreach($entry in $controllers.ToArray()){
            $entry.Process.Refresh()
            $row=$state.models[$entry.Model-1]
            if($entry.Process.HasExited){$entry.Process.WaitForExit();$entry.Process.Refresh();$row.status=if($entry.Process.ExitCode -eq 0){'completed'}else{'failed'};if($entry.Process.ExitCode -ne 0){$failed+=$entry.Model}}
            else{$active++;$row.status='running'}
        }
        $state.updated_utc=[DateTime]::UtcNow.ToString('o');Save-State $statePath $state
        if(-not $active){if($failed.Count){throw "Failed model controllers: $($failed -join ', ')"};break}
        Start-Sleep -Seconds 5
    }
    if([int]$control[5] -eq 1){
        $merged=Join-Path $run 'merged'
        New-Item -ItemType Directory -Path $merged -Force|Out-Null
        $allDist=@();$allTargets=@()
        foreach($model in 1..5){
            $result=Join-Path $run ("model-$model/results/$name")
            $resultState=Read-State (Join-Path $result 'manifest.json')
            if($resultState.stages.A -ne 'completed' -or $resultState.stages.merge -ne 'completed'){throw "Model $model analysis is incomplete."}
            $dest=Join-Path $merged ('model-'+$model)
            if(-not(Test-Path -LiteralPath $dest)){Copy-Item -LiteralPath (Join-Path $result "merged/model-$model") -Destination $dest -Recurse}
            $allDist+=@(Import-Csv -LiteralPath (Join-Path $result 'merged/distributions.csv'))
            $allTargets+=@(Import-Csv -LiteralPath (Join-Path $result 'merged/target_isotopes.csv'))
        }
        foreach($pair in @(@('distributions.csv',$allDist),@('target_isotopes.csv',$allTargets))){
            $dest=Join-Path $merged $pair[0]
            if(-not(Test-Path -LiteralPath $dest)){[IO.File]::WriteAllLines($dest,[string[]]@($pair[1]|ConvertTo-Csv -NoTypeInformation),$utf8)}
        }
    }
    $state.stage='completed';$state.updated_utc=[DateTime]::UtcNow.ToString('o');Save-State $statePath $state
    Write-Output "Completed five table batch: $run"
}catch{
    if($state -and $statePath){$state.stage='failed';$state['error']=$_.Exception.Message;$state.updated_utc=[DateTime]::UtcNow.ToString('o');Save-State $statePath $state}
    throw
}finally{
    if($locked){$mutex.ReleaseMutex()};$mutex.Dispose()
}
