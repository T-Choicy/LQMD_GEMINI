$ErrorActionPreference='Stop'
$repo=Split-Path -Parent $PSScriptRoot
$run=Join-Path $repo ('tmp/tests/example/'+[guid]::NewGuid().ToString('N'))
$paths=@(Get-ChildItem (Join-Path $repo 'fortran/Input file') -File)
$before=@{};foreach($file in $paths){$before[$file.Name]=(Get-FileHash $file.FullName).Hash}
& (Join-Path $repo 'scripts/run_example.ps1') -PrepareOnly -EnergyMeV 40 -Shards 2 -EventsPerShard 3 -Model 2 -RunDirectory $run
if($LASTEXITCODE -ne 0){throw 'Preparation returned a nonzero status.'}
$m=Get-Content (Join-Path $run 'run.json') -Raw|ConvertFrom-Json
if($m.energy_MeV -ne 40 -or $m.total_events -ne 6 -or $m.model -ne 2 -or $m.exit_code -ne 0){throw 'Run metadata mismatch.'}
$energy=([IO.File]::ReadAllLines((Join-Path $run 'source/Input file/LQMD0.in'))[3] -split ',')[0]
if([double]::Parse($energy,[Globalization.CultureInfo]::InvariantCulture) -ne 0.04){throw 'Beam-energy conversion mismatch.'}
try{& (Join-Path $repo 'scripts/run_example.ps1') -RunDirectory $run;throw 'Existing directory accepted.'}catch{if($_.Exception.Message -notlike '*RunDirectory already exists*'){throw}}
try{& (Join-Path $repo 'scripts/run_example.ps1') -Shards 0;throw 'Invalid shards accepted.'}catch{if($_.Exception.Message -notmatch 'Shards'){throw}}
foreach($file in $paths){if((Get-FileHash $file.FullName).Hash -ne $before[$file.Name]){throw "Shared input modified: $($file.Name)"}}
Write-Output 'PASS example preparation, energy/metadata, event partition, existing-directory rejection, argument validation and source preservation.'
