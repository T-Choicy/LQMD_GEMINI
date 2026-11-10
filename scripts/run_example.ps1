[CmdletBinding()]
param(
    [ValidateRange(1,25)][int]$Shards=1,
    [ValidateRange(1,2147483647)][int]$EventsPerShard=2,
    [ValidateRange(1,5)][int]$Model=1,
    [ValidateRange(1,100000)][double]$EnergyMeV=80,
    [string]$RunDirectory,
    [switch]$PrepareOnly
)
$ErrorActionPreference='Stop'
$repo=Split-Path -Parent $PSScriptRoot
if([long]$Shards*$EventsPerShard -gt [int]::MaxValue){throw 'Total events exceed the supported range.'}
$energy=$EnergyMeV.ToString('0.########',[Globalization.CultureInfo]::InvariantCulture)
if(-not $RunDirectory){$RunDirectory=Join-Path $repo ('runs/p232Th'+$energy+'-model'+$Model+'/'+[DateTime]::UtcNow.ToString('yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8))}
$run=[IO.Path]::GetFullPath($RunDirectory)
if(Test-Path -LiteralPath $run){throw 'RunDirectory already exists; use a new directory to preserve previous results.'}
$binary=Join-Path $repo 'fortran/x64/Release/win_batch'
foreach($name in @('ctrlL_win.exe','LQMD.exe','exfrag.exe')){if(-not(Test-Path (Join-Path $binary $name))){throw 'Build the programs first: ./scripts/build_windows_controllers.ps1'}}
$source=Join-Path $run 'source'
New-Item -ItemType Directory $source -Force|Out-Null
foreach($name in @('Input file','gdat','tbl','pini.dat','tini.dat')){Copy-Item -LiteralPath (Join-Path $repo ('fortran/'+$name)) -Destination $source -Recurse}
foreach($name in @('LQMD0.in','anahyp0.in','CNPCS_in.in')){Copy-Item -LiteralPath (Join-Path $repo ('examples/p232Th-80MeV/'+$name)) -Destination (Join-Path $source ('Input file/'+$name)) -Force}
$utf8=[Text.UTF8Encoding]::new($false)
$physics=Join-Path $source 'Input file/LQMD0.in'
$lines=[IO.File]::ReadAllLines($physics)
$gev=($EnergyMeV/1000).ToString('0.########',[Globalization.CultureInfo]::InvariantCulture)
$lines[3]=[regex]::Replace($lines[3],'^\s*[+-]?[0-9.]+(?:[EeDd][+-]?\d+)?',$gev)
[IO.File]::WriteAllLines($physics,$lines,$utf8)
[IO.File]::WriteAllLines((Join-Path $source 'Input file/INPUT.in'),[string[]]@("$Shards","$EventsPerShard",'calculation','0','1','1','0'),$utf8)
[IO.File]::WriteAllLines((Join-Path $source 'Input file/m_select.in'),[string[]]@("$Model",'1','314159'),$utf8)
$commit='unavailable';$dirty=$null
if(Get-Command git -ErrorAction SilentlyContinue){
 $value=& git -C $repo rev-parse HEAD 2>$null
 if($LASTEXITCODE -eq 0){$commit=([string]$value).Trim();$dirty=[bool](& git -C $repo status --porcelain)}
}
$hashes=[ordered]@{}
foreach($name in @('ctrlL_win.exe','LQMD.exe','exfrag.exe')){$hashes[$name]=(Get-FileHash (Join-Path $binary $name)).Hash}
$record=[ordered]@{version=([IO.File]::ReadAllText((Join-Path $repo 'VERSION'))).Trim();commit=$commit;dirty_worktree=$dirty;started_utc=[DateTime]::UtcNow.ToString('o');projectile='p';target='232Th';energy_MeV=$EnergyMeV;model=$Model;mass_framework=1;shards=$Shards;events_per_shard=$EventsPerShard;total_events=([long]$Shards*$EventsPerShard);base_seed=314159;prepare_only=[bool]$PrepareOnly;executable_sha256=$hashes}
$meta=Join-Path $run 'run.json'
[IO.File]::WriteAllText($meta,($record|ConvertTo-Json -Depth 6),$utf8)
$params=@{FilePath=(Join-Path $binary 'ctrlL_win.exe');WorkingDirectory=$source;WindowStyle='Hidden';PassThru=$true;RedirectStandardOutput=(Join-Path $run 'controller.stdout.log');RedirectStandardError=(Join-Path $run 'controller.stderr.log')}
if($PrepareOnly){$params.ArgumentList='--prepare-only'}
$process=Start-Process @params
$null=$process.Handle
$process.WaitForExit();$process.Refresh()
$record.exit_code=$process.ExitCode;$record.finished_utc=[DateTime]::UtcNow.ToString('o')
[IO.File]::WriteAllText($meta,($record|ConvertTo-Json -Depth 6),$utf8)
if($process.ExitCode -ne 0){throw "Calculation failed ($($process.ExitCode)); inspect $run/controller.stderr.log and shard logs."}
if($PrepareOnly){Write-Output "Prepared inputs: $source/results/calculation"}else{Write-Output "Results: $source/results/calculation/merged"}
Write-Output "Run provenance: $meta"
exit 0
