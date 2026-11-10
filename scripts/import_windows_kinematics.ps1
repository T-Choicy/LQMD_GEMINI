[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$PreviousRunDirectory,
    [Parameter(Mandatory=$true)][string]$RunDirectory
)
$ErrorActionPreference='Stop'
$env:PSModulePath=(Join-Path $PSHOME 'Modules')+[IO.Path]::PathSeparator+$env:PSModulePath
Import-Module Microsoft.PowerShell.Management
Import-Module Microsoft.PowerShell.Utility
$old=[IO.Path]::GetFullPath($PreviousRunDirectory).TrimEnd('\')
$new=[IO.Path]::GetFullPath($RunDirectory).TrimEnd('\')
if($old -eq $new){throw 'Import requires a separate prepared batch.'}
$utf8=New-Object Text.UTF8Encoding($false)
function Read-Json([string]$path){return ([IO.File]::ReadAllText($path)|ConvertFrom-Json)}
function Save-Json([string]$path,$value){[IO.File]::WriteAllText($path,($value|ConvertTo-Json -Depth 15),$utf8)}
function Same-Map($left,$right){
    $names=@($left.PSObject.Properties.Name)
    if($names.Count -ne @($right.PSObject.Properties.Name).Count){return $false}
    foreach($name in $names){if([string]$left.$name -cne [string]$right.$name){return $false}}
    return $true
}
$sourceBatch=Read-Json (Join-Path $old 'batch-manifest.json')
$targetBatch=Read-Json (Join-Path $new 'batch-manifest.json')
foreach($key in @('events_per_model','tasks_per_model')){
    if($sourceBatch.identity.$key -ne $targetBatch.identity.$key){throw "Import layout differs: $key"}
}
$rows=New-Object 'System.Collections.Generic.List[object]'
foreach($model in 1..5){
    $sourceRoot=Join-Path $old "model-$model"
    $targetRoot=Join-Path $new "model-$model"
    $control=Join-Path $sourceRoot 'Input file/INPUT.in'
    $controlRows=@([IO.File]::ReadAllLines($control)|ForEach-Object{($_ -split '[!#]',2)[0].Trim()}|Where-Object{$_})
    $name=$controlRows[2]
    $sourceRun=Join-Path $sourceRoot ('results/'+$name)
    $targetRun=Join-Path $targetRoot ('results/'+$name)
    if(Test-Path -LiteralPath (Join-Path $targetRun 'manifest.json')){throw "Import target already initialized: $targetRun"}
    $sourceManifestPath=Join-Path $sourceRun 'manifest.json'
    $sourceManifest=Read-Json $sourceManifestPath
    if($sourceManifest.stages.L -ne 'completed'){throw "Source dynamics incomplete: model $model"}
    $prepare=Start-Process -FilePath (Join-Path $targetRoot 'ctrlL_win.exe') -ArgumentList '--prepare-only' -WorkingDirectory $targetRoot -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $targetRoot 'import-prepare.stdout.log') -RedirectStandardError (Join-Path $targetRoot 'import-prepare.stderr.log')
    $null=$prepare.Handle;$prepare.WaitForExit();$prepare.Refresh()
    if($prepare.ExitCode -ne 0){throw "Import preparation failed: model $model"}
    $targetManifestPath=Join-Path $targetRun 'manifest.json'
    $targetManifest=Read-Json $targetManifestPath
    if(($sourceManifest.shard_plan|ConvertTo-Json -Compress) -cne ($targetManifest.shard_plan|ConvertTo-Json -Compress)){throw "Event ranges or seeds differ: model $model"}
    foreach($shard in $sourceManifest.shard_plan){
        $label='shard-{0:d4}' -f [int]$shard.Index
        $valid=@()
        foreach($attempt in @(Get-ChildItem -LiteralPath (Join-Path $sourceRun "kinematics/$label") -Directory -Filter 'attempt-*'|Sort-Object Name -Descending)){
            $sourceState=Read-Json (Join-Path $attempt.FullName 'state.json')
            if($sourceState.status -eq 'completed'){$valid=@($attempt);break}
        }
        if(-not $valid.Count){throw "No completed source: model $model $label"}
        $sourceAttempt=$valid[0].FullName
        $targetAttempt=Join-Path $targetRun "kinematics/$label/attempt-0001"
        $sourceStatePath=Join-Path $sourceAttempt 'state.json'
        $sourceState=Read-Json $sourceStatePath
        $targetState=Read-Json (Join-Path $targetAttempt 'state.json')
        if($sourceState.exit_code -ne 0 -or $sourceState.events -ne $shard.Events){throw 'Invalid source completion metadata.'}
        if($sourceState.first_event -ne $shard.First -or $sourceState.last_event -ne $shard.Last -or $sourceState.seed -ne $shard.Seed){throw 'Source event range or seed differs from the shard plan.'}
        if(-not(Same-Map $sourceState.input_hashes $targetState.input_hashes)){throw "Physical inputs differ: model $model $label"}
        $sourceTask=Read-Json (Join-Path $sourceAttempt 'task.json')
        $targetTask=Read-Json (Join-Path $targetAttempt 'task.json')
        if($sourceTask.ExecutableHash -cne $targetTask.ExecutableHash){throw 'Dynamical executable changed.'}
        foreach($entry in $sourceState.input_hashes.PSObject.Properties){
            if((Get-FileHash -LiteralPath (Join-Path $sourceAttempt $entry.Name) -Algorithm SHA256).Hash -cne $entry.Value){throw "Source input changed: $($entry.Name)"}
        }
        foreach($file in @('coll.dat','initangm.dat')){
            $path=Join-Path $sourceAttempt $file
            if((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -cne [string]$sourceState.output_hashes.$file){throw "Source output changed: $path"}
            # Both trees are immutable source data; analysis only reads these links.
            New-Item -ItemType HardLink -Path (Join-Path $targetAttempt $file) -Target $path|Out-Null
        }
        foreach($file in @('program.stdout.log','program.stderr.log')){Copy-Item -LiteralPath (Join-Path $sourceAttempt $file) -Destination $targetAttempt}
        Copy-Item -LiteralPath $sourceStatePath -Destination (Join-Path $targetAttempt 'imported-source-state.json')
        Copy-Item -LiteralPath (Join-Path $sourceAttempt 'task.json') -Destination (Join-Path $targetAttempt 'imported-source-task.json')
        $lineage=[ordered]@{source_attempt=$sourceAttempt;source_state_sha256=(Get-FileHash -LiteralPath $sourceStatePath).Hash;source_manifest_sha256=(Get-FileHash -LiteralPath $sourceManifestPath).Hash;imported_utc=[DateTime]::UtcNow.ToString('o');method='validated NTFS hard links; no LQMD execution'}
        $sourceState|Add-Member -NotePropertyName imported_kinematics -NotePropertyValue $lineage -Force
        Save-Json (Join-Path $targetAttempt 'state.json') $sourceState
        $rows.Add([ordered]@{model=$model;shard=$shard.Index;source=$sourceAttempt;destination=$targetAttempt;output_hashes=$sourceState.output_hashes})
    }
    $targetManifest.stages.L='completed'
    $targetManifest|Add-Member -NotePropertyName kinematics_import -NotePropertyValue ([ordered]@{source_run=$sourceRun;source_manifest_sha256=(Get-FileHash -LiteralPath $sourceManifestPath).Hash;imported_utc=[DateTime]::UtcNow.ToString('o');notice='Dynamical events were computed in the source run, not recomputed here.'}) -Force
    Save-Json $targetManifestPath $targetManifest
    Write-Output "Imported validated model $model dynamics; analysis remains pending."
}
Save-Json (Join-Path $new 'kinematics-import.json') ([ordered]@{previous_run=$old;derived_run=$new;shards=$rows.ToArray();completed_utc=[DateTime]::UtcNow.ToString('o')})
