[CmdletBinding()]
param()
$ErrorActionPreference='Stop'
$env:PSModulePath=(Join-Path $PSHOME 'Modules')+[IO.Path]::PathSeparator+$env:PSModulePath
Import-Module Microsoft.PowerShell.Utility
Import-Module Microsoft.PowerShell.Management
$repo=Split-Path -Parent $PSScriptRoot
$test=Join-Path $repo ('tmp/tests/windows-legacy/'+[guid]::NewGuid().ToString('N'))
$source=Join-Path $test 'source'
$utf8=New-Object Text.UTF8Encoding($false)
New-Item -ItemType Directory -Path $source -Force|Out-Null
foreach($folder in @('Input file','gdat','tbl')){Copy-Item -LiteralPath (Join-Path $repo ('fortran/'+$folder)) -Destination $source -Recurse}
foreach($file in @('pini.dat','tini.dat')){Copy-Item -LiteralPath (Join-Path $repo ('fortran/'+$file)) -Destination $source}
[IO.File]::WriteAllLines((Join-Path $source 'Input file/INPUT.in'),[string[]]@('5 !NFILE','2 !NRUN_PER','mock_25 !NAME','0 !T_1ST','1 !T_PER','1 !CALL_ANA','0 !FAIL'),$utf8)
$binary=Join-Path $source 'x64/Release/win_batch'
New-Item -ItemType Directory -Path $binary -Force|Out-Null
foreach($file in @('ctrlL_win.exe','ctrlA_win.exe')){Copy-Item -LiteralPath (Join-Path $repo ('fortran/x64/Release/win_batch/'+$file)) -Destination $binary}
$compiler='C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\ifx.exe'
$setvars='C:\Program Files (x86)\Intel\oneAPI\setvars.bat'
$fixture=Join-Path $repo 'tests/fixtures/windows_controller_mock.f90'
$build=Join-Path $test 'mock-build.cmd'
[IO.File]::WriteAllLines($build,[string[]]@('@echo off',('call "'+$setvars+'" >nul'),('"'+$compiler+'" /nologo /free /O2 /libs:static /object:"'+(Join-Path $test 'mock.obj')+'" /exe:"'+(Join-Path $binary 'LQMD.exe')+'" "'+$fixture+'"'),'exit /b %errorlevel%'),[Text.Encoding]::ASCII)
$buildProc=Start-Process -FilePath $env:ComSpec -ArgumentList ('/d /s /c ""'+$build+'""') -WorkingDirectory $test -PassThru -WindowStyle Hidden -RedirectStandardOutput (Join-Path $test 'build.stdout.log') -RedirectStandardError (Join-Path $test 'build.stderr.log')
$null=$buildProc.Handle
$buildProc.WaitForExit();$buildProc.Refresh()
if($buildProc.ExitCode -ne 0){throw ('ifx mock compilation failed: '+[IO.File]::ReadAllText((Join-Path $test 'build.stderr.log'))+' '+[IO.File]::ReadAllText((Join-Path $test 'build.stdout.log')))}
Copy-Item -LiteralPath (Join-Path $binary 'LQMD.exe') -Destination (Join-Path $binary 'exfrag.exe')
$run=Join-Path $test 'five model run'
$release=Join-Path $test 'release.mock'
$beforeAllow=$env:MOCK_ALLOW_SELECTED_LQMD
$beforeRelease=$env:MOCK_RELEASE_FILE
$env:MOCK_ALLOW_SELECTED_LQMD='1';$env:MOCK_RELEASE_FILE=$release
$controller=$null
try{
    $launcher=Join-Path $repo 'scripts/run_windows_mass_models.ps1'
    $args='-NoProfile -ExecutionPolicy Bypass -File "'+$launcher+'" -SourceDirectory "'+$source+'" -RunDirectory "'+$run+'"'
    $controller=Start-Process -FilePath powershell.exe -ArgumentList $args -PassThru -WindowStyle Hidden -RedirectStandardOutput (Join-Path $test 'batch.stdout.log') -RedirectStandardError (Join-Path $test 'batch.stderr.log')
    $null=$controller.Handle
    $deadline=[DateTime]::UtcNow.AddSeconds(55)
    $count=0
    do{
        $count=0
        $perModel=@{}
        foreach($process in @(Get-Process -Name LQMD -ErrorAction SilentlyContinue)){
            if($process.Path -and $process.Path.StartsWith($run,[StringComparison]::OrdinalIgnoreCase)){
                $count++
                $parent=Split-Path -Leaf (Split-Path -Parent $process.Path)
                if(-not $perModel.ContainsKey($parent)){$perModel[$parent]=0};$perModel[$parent]++
            }
        }
        $controller.Refresh()
        if($controller.HasExited){throw ('Batch ended before 25 native processes: '+[IO.File]::ReadAllText((Join-Path $test 'batch.stderr.log')))}
        if($count -eq 25){break}
        Start-Sleep -Milliseconds 150
    }while([DateTime]::UtcNow -lt $deadline)
    if($count -ne 25){throw "Expected 25 simultaneous LQMD processes, found $count. See $test"}
    foreach($model in 1..5){if($perModel['model-'+$model] -ne 5){throw "Expected 5 LQMD for model $model"}}
    [IO.File]::WriteAllText($release,'release',$utf8)
    $controller.WaitForExit();$controller.Refresh()
    if($controller.ExitCode -ne 0){throw ('Batch failed: '+[IO.File]::ReadAllText((Join-Path $test 'batch.stderr.log')))}
    foreach($model in 1..5){
        $modelRoot=Join-Path $run ('model-'+$model)
        $result=Join-Path $modelRoot 'results/mock_25'
        $manifest=Get-Content -LiteralPath (Join-Path $result 'manifest.json') -Raw|ConvertFrom-Json
        if($manifest.events -ne 10 -or $manifest.shards -ne 5 -or $manifest.stages.L -ne 'completed' -or $manifest.stages.A -ne 'completed'){throw "Wrong event counts/auto-analysis model $model"}
        $ids=@()
        foreach($shard in 1..5){
            $piece=Join-Path $result ('kinematics/shard-{0:d4}/attempt-0001' -f $shard)
            foreach($line in [IO.File]::ReadAllLines((Join-Path $piece 'initangm.dat'))){$fields=$line.Trim() -split '\s+';$ids+=[int]$fields[0];if([int]$fields[2] -ne $model){throw 'LQMD selected model mismatch'}}
        }
        if(($ids|Sort-Object) -join ',' -ne '1,2,3,4,5,6,7,8,9,10'){throw 'Event IDs repeated or missing'}
        foreach($name in @('LQMD.in','anahyp.in')){
            if((Get-FileHash -LiteralPath (Join-Path $source ('Input file/'+$name))).Hash -ne (Get-FileHash -LiteralPath (Join-Path $modelRoot ('Input file/'+$name))).Hash){throw 'Original actual input was changed'}
        }
    }
    if(@(Import-Csv -LiteralPath (Join-Path $run 'merged/distributions.csv')).Count -ne 6000 -or @(Import-Csv -LiteralPath (Join-Path $run 'merged/target_isotopes.csv')).Count -ne 25){throw 'Combined CSV row counts invalid'}
    Write-Output "PASS: 25 simultaneous native mocks, five per table; original inputs, exact per-table events, model-isolated data and automatic analysis: $test"
}finally{
    if(-not(Test-Path -LiteralPath $release)){[IO.File]::WriteAllText($release,'release',$utf8)}
    if($controller -and -not $controller.HasExited){$controller.WaitForExit(30000)|Out-Null}
    if($null -eq $beforeAllow){Remove-Item Env:MOCK_ALLOW_SELECTED_LQMD -ErrorAction SilentlyContinue}else{$env:MOCK_ALLOW_SELECTED_LQMD=$beforeAllow}
    if($null -eq $beforeRelease){Remove-Item Env:MOCK_RELEASE_FILE -ErrorAction SilentlyContinue}else{$env:MOCK_RELEASE_FILE=$beforeRelease}
}
