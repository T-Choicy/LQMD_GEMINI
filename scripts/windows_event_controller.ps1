[CmdletBinding()]
param(
    [ValidateSet('L','A')][string]$Stage,
    [string]$Config,
    [switch]$PrepareOnly,
    [switch]$MergeOnly,
    [string]$InternalTask
)

$ErrorActionPreference = 'Stop'
$script:WindowsHost = [Environment]::OSVersion.Platform -eq [PlatformID]::Win32NT
$script:ChildPowerShell = if ($script:WindowsHost) { 'powershell.exe' } else { Join-Path $PSHOME 'pwsh' }
$script:HeldSlots = @{}
$requiredModules = Join-Path $PSHOME 'Modules'
if (($env:PSModulePath -split [IO.Path]::PathSeparator) -notcontains $requiredModules) {
    $env:PSModulePath = $requiredModules + [IO.Path]::PathSeparator + $env:PSModulePath
}
Import-Module Microsoft.PowerShell.Management -ErrorAction Stop
Import-Module Microsoft.PowerShell.Utility -ErrorAction Stop
$script:Utf8NoBom = New-Object System.Text.UTF8Encoding($false)
$script:Models = @{
    '1' = 'FRDM12'; '2' = 'HFB02'; '3' = 'KTUY05'; '4' = 'WS4+RBF'; '5' = 'LDMM'
}
$script:targets = @(
    @{ A = 225; Z = 89; N = 136; Name = '225Ac' },
    @{ A = 227; Z = 89; N = 138; Name = '227Ac' },
    @{ A = 227; Z = 90; N = 137; Name = '227Th' },
    @{ A = 225; Z = 88; N = 137; Name = '225Ra' },
    @{ A = 223; Z = 88; N = 135; Name = '223Ra' }
)
. (Join-Path $PSScriptRoot 'windows_controller_results.ps1')

function Write-JsonAtomic {
    param([string]$Path, [object]$Value)
    $tmp = $Path + '.tmp-' + [guid]::NewGuid().ToString('N')
    $json = $Value | ConvertTo-Json -Depth 20
    [IO.File]::WriteAllText($tmp, $json + "`n", $script:Utf8NoBom)
    if (Test-Path -LiteralPath $Path) { throw "Refusing to overwrite existing state file: $Path" }
    [IO.File]::Move($tmp, $Path)
}

function Set-JsonAtomic {
    param([string]$Path, [object]$Value)
    $tmp = $Path + '.tmp-' + [guid]::NewGuid().ToString('N')
    $json = $Value | ConvertTo-Json -Depth 20
    [IO.File]::WriteAllText($tmp, $json + "`n", $script:Utf8NoBom)
    for ($attempt=0; ; $attempt++) {
        $backup = $Path + '.bak-' + [guid]::NewGuid().ToString('N')
        try {
            if ([IO.File]::Exists($Path)) { [IO.File]::Replace($tmp, $Path, $backup) }
            else { [IO.File]::Move($tmp, $Path) }
            return
        } catch [IO.IOException] {
            if ($attempt -ge 39) { throw }
            Start-Sleep -Milliseconds 50
        } finally { if ([IO.File]::Exists($backup)) { [IO.File]::Delete($backup) } }
    }
}

function Read-JsonShared {
    param([string]$Path)
    $stream = [IO.File]::Open($Path,[IO.FileMode]::Open,[IO.FileAccess]::Read,([IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete))
    $reader = New-Object IO.StreamReader($stream, [Text.Encoding]::UTF8, $true)
    try { return ($reader.ReadToEnd() | ConvertFrom-Json) }
    finally { $reader.Dispose(); $stream.Dispose() }
}

function Get-Sha256 {
    param([string]$Path)
    $stream = [IO.File]::Open($Path,[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete)
    $sha = [Security.Cryptography.SHA256]::Create()
    try { return ([BitConverter]::ToString($sha.ComputeHash($stream))).Replace('-','').ToUpperInvariant() }
    finally { $stream.Dispose(); $sha.Dispose() }
}

function Resolve-ConfigPath {
    param([string]$Base, [string]$Value)
    if ([IO.Path]::IsPathRooted($Value)) { return [IO.Path]::GetFullPath($Value) }
    return [IO.Path]::GetFullPath((Join-Path $Base $Value))
}

function Get-Config {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { throw "Configuration not found: $Path" }
    $configPath = [IO.Path]::GetFullPath($Path)
    $base = Split-Path -Parent $configPath
    if ([IO.Path]::GetExtension($configPath).Equals('.in',[StringComparison]::OrdinalIgnoreCase)) {
        $leaf = [IO.Path]::GetFileName($configPath)
        if ($leaf -notin @('INPUT.in','INPUT_ANA.in')) { throw 'Legacy controller input must be named INPUT.in or INPUT_ANA.in.' }
        $source = Split-Path -Parent $base
        $inputRoot = Join-Path $source 'Input file'
        function Read-LegacyRows([string]$file, [int]$count) {
            $rows = New-Object 'System.Collections.Generic.List[string]'
            foreach ($raw in [IO.File]::ReadAllLines($file)) {
                $line = ($raw -split '[!#]',2)[0].Trim()
                if ($line) { $rows.Add($line); if ($rows.Count -ge $count) { break } }
            }
            if ($rows.Count -lt $count) { throw "Expected $count data rows in $file." }
            return ,$rows.ToArray()
        }
        $isAnalysis = $leaf -eq 'INPUT_ANA.in'
        if ($isAnalysis) { $rows = Read-LegacyRows $configPath 9; $nfile=[int]$rows[0]; $per=[int]$rows[1]; $name=$rows[2]; $tper=$rows[3]; $xcut=$rows[4] }
        else { $rows = Read-LegacyRows $configPath 7; $nfile=[int]$rows[0]; $per=[int]$rows[1]; $name=$rows[2]; $tper=$rows[4] }
        if ($nfile -lt 1 -or $nfile -gt 25 -or $per -lt 1) { throw 'Legacy NFILE must be 1..25 and NRUN_PER must be positive.' }
        if ($name -in @('.','..') -or $name.IndexOfAny([IO.Path]::GetInvalidFileNameChars()) -ge 0 -or $name.Contains('/') -or $name.Contains('\')) { throw 'NAME must be a plain task name, without directory separators.' }
        if ([long]$nfile*$per -gt [int]::MaxValue) { throw 'NFILE * NRUN_PER exceeds the program event range.' }
        if (-not $isAnalysis -and [int]$rows[6] -ne 0) { throw 'Legacy FAIL must be 0; failed-event offsets are not supported.' }
        $selectorPath = Join-Path $inputRoot 'm_select.in'
        $selRows = Read-LegacyRows $selectorPath 3
        $model = [int]$selRows[0]; $framework = [int]$selRows[1]; $seed=[int]$selRows[2]
        if ($model -lt 1 -or $model -gt 5) { throw 'Legacy m_select model must be 1..5.' }
        if ($framework -ne 1) { throw 'Legacy runs require mass_framework=1 to keep LQMD dynamics frozen at M0.' }
        $repo = Split-Path -Parent $PSScriptRoot
        $exeRoot = Join-Path $repo 'fortran/x64/Release/win_batch'
        if (-not $script:WindowsHost) { $exeRoot = Join-Path $repo 'fortran/bin' }
        $cfg = [pscustomobject]@{ SourceDirectory=$source; RunDirectory=(Join-Path $source ('results/' + $name)); LqmdExecutable=(Join-Path $exeRoot 'LQMD.exe'); ExfragExecutable=(Join-Path $exeRoot 'exfrag.exe'); Events=($nfile*$per); Shards=$nfile; MaxWorkers=$nfile; BaseSeed=$seed; Models=@($model); SelectedModel=$model; Legacy=$true; LegacyAnalysis=$isAnalysis; ConfigPath=$configPath; ConfigBase=$base; LegacyRows=$rows }
        foreach ($field in @('LqmdExecutable','ExfragExecutable')) {
            if (-not $script:WindowsHost) { $cfg.$field = $cfg.$field -replace '\.exe$','' }
            $local = Join-Path $source ([IO.Path]::GetFileName($cfg.$field))
            if (Test-Path -LiteralPath $local -PathType Leaf) { $cfg.$field=$local }
        }
        foreach ($field in @('LqmdExecutable','ExfragExecutable')) { if (-not (Test-Path -LiteralPath $cfg.$field -PathType Leaf)) { throw "Executable not found: $($cfg.$field)" } }
        foreach ($required in @('Input file','gdat','tbl','pini.dat','tini.dat')) { if (-not (Test-Path -LiteralPath (Join-Path $source $required))) { throw "Required input missing from SourceDirectory: $required" } }
        if ($isAnalysis) {
            $cfg | Add-Member -NotePropertyName LegacyAnalysisRows -NotePropertyValue $rows -Force
            $cfg | Add-Member -NotePropertyName LegacyAnalysisXcut -NotePropertyValue $xcut -Force
        } else {
            if ([int]$rows[5] -notin @(0,1)) { throw 'CALL_ANA must be 0 or 1.' }
            $cfg | Add-Member -NotePropertyName CallAnalysis -NotePropertyValue ([int]$rows[5] -eq 1) -Force
            if ($cfg.CallAnalysis) { Sync-LegacyAnalysisControl $cfg }
        }
        return $cfg
    }
    $cfg = Get-Content -LiteralPath $configPath -Raw -Encoding UTF8 | ConvertFrom-Json
    foreach ($name in @('SourceDirectory','RunDirectory','LqmdExecutable','ExfragExecutable')) {
        if (-not $cfg.PSObject.Properties[$name] -or [string]::IsNullOrWhiteSpace([string]$cfg.$name)) {
            throw "Missing required JSON field: $name"
        }
    }
    $cfg | Add-Member -NotePropertyName ConfigPath -NotePropertyValue $configPath -Force
    $cfg | Add-Member -NotePropertyName ConfigBase -NotePropertyValue $base -Force
    foreach ($name in @('SourceDirectory','RunDirectory','LqmdExecutable','ExfragExecutable')) {
        $cfg.$name = Resolve-ConfigPath $base ([string]$cfg.$name)
    }
    foreach ($name in @('Events','Shards','MaxWorkers','BaseSeed','Models')) {
        if (-not $cfg.PSObject.Properties[$name]) { throw "Missing required JSON field: $name" }
    }
    $cfg.Events = [int]$cfg.Events
    $cfg.Shards = [int]$cfg.Shards
    $cfg.MaxWorkers = [int]$cfg.MaxWorkers
    $cfg.BaseSeed = [int]$cfg.BaseSeed
    $cfg.Models = @($cfg.Models | ForEach-Object { [int]$_ })
    if ($cfg.Events -lt 1 -or $cfg.Shards -lt 1 -or $cfg.Shards -gt $cfg.Events) { throw 'Events and Shards must be positive, with Shards <= Events.' }
    if ($cfg.MaxWorkers -lt 1 -or $cfg.MaxWorkers -gt 8) { throw 'MaxWorkers must be in the safe range 1..8.' }
    if ($cfg.BaseSeed -lt 1 -or $cfg.BaseSeed -gt 2000000000) { throw 'BaseSeed must be in the range 1..2000000000.' }
    if ($cfg.Models.Count -ne 5 -or (($cfg.Models | Sort-Object -Unique) -join ',') -ne '1,2,3,4,5') { throw 'Models must be exactly [1,2,3,4,5].' }
    if (-not (Test-Path -LiteralPath $cfg.SourceDirectory -PathType Container)) { throw "SourceDirectory not found: $($cfg.SourceDirectory)" }
    foreach ($name in @('LqmdExecutable','ExfragExecutable')) { if (-not (Test-Path -LiteralPath $cfg.$name -PathType Leaf)) { throw "Executable not found: $($cfg.$name)" } }
    $src = $cfg.SourceDirectory
    foreach ($name in @('Input file','gdat','tbl','pini.dat','tini.dat')) {
        if (-not (Test-Path -LiteralPath (Join-Path $src $name))) { throw "Required input missing from SourceDirectory: $name" }
    }
    foreach ($folder in @('Input file','gdat','tbl')) {
        $protected = [IO.Path]::GetFullPath((Join-Path $src $folder)).TrimEnd([IO.Path]::DirectorySeparatorChar)
        $comparison = if ($script:WindowsHost) {[StringComparison]::OrdinalIgnoreCase} else {[StringComparison]::Ordinal}
        if ($cfg.RunDirectory.Equals($protected,$comparison) -or $cfg.RunDirectory.StartsWith($protected+[IO.Path]::DirectorySeparatorChar,$comparison)) { throw 'RunDirectory must be outside the source input/table directories.' }
    }
    return $cfg
}

function Get-InputSnapshot {
    param([string]$Root, [string[]]$ExecutablePaths, [switch]$SkipInitialFiles)
    $map = [ordered]@{}
    $files = New-Object 'System.Collections.Generic.List[System.IO.FileInfo]'
    foreach ($folder in @('Input file','gdat','tbl')) {
        Get-ChildItem -LiteralPath (Join-Path $Root $folder) -File -Recurse | ForEach-Object { $files.Add($_) }
    }
    if (-not $SkipInitialFiles) { foreach ($name in @('pini.dat','tini.dat')) { $files.Add((Get-Item -LiteralPath (Join-Path $Root $name))) } }
    $rootPrefix = [IO.Path]::GetFullPath($Root).TrimEnd([IO.Path]::DirectorySeparatorChar) + [IO.Path]::DirectorySeparatorChar
    foreach ($file in $files | Sort-Object FullName) {
        $rel = $file.FullName.Substring($rootPrefix.Length).Replace([string][IO.Path]::DirectorySeparatorChar, '/')
        $map[$rel] = Get-Sha256 $file.FullName
    }
    foreach ($exe in $ExecutablePaths) { $map['exe/' + [IO.Path]::GetFileName($exe)] = Get-Sha256 $exe }
    return $map
}

function Copy-RunInputs {
    param([string]$Source, [string]$Destination, [switch]$UseTemplates)
    New-Item -ItemType Directory -Path $Destination -Force | Out-Null
    foreach ($folder in @('Input file','gdat','tbl')) {
        Copy-Item -LiteralPath (Join-Path $Source $folder) -Destination $Destination -Recurse
    }
    foreach ($name in @('pini.dat','tini.dat')) { Copy-Item -LiteralPath (Join-Path $Source $name) -Destination $Destination }
    if ($UseTemplates) { foreach ($name in @('LQMD.in','anahyp.in')) {
        $template = Join-Path $Source ('Input file/' + [IO.Path]::GetFileNameWithoutExtension($name) + '0.in')
        if (-not (Test-Path -LiteralPath $template -PathType Leaf)) { throw "Required physical template not found: $template" }
        Copy-Item -LiteralPath $template -Destination (Join-Path $Destination ('Input file/' + $name)) -Force
    } }
}

function Sync-LegacyAnalysisControl {
    param($Cfg)
    $rows = @([IO.File]::ReadAllLines((Join-Path $Cfg.SourceDirectory 'Input file/LQMD0.in')) | ForEach-Object {($_ -split '!',2)[0].Trim()} | Where-Object { $_ -and -not $_.StartsWith('#') })
    if ($rows.Count -lt 13) { throw 'LQMD0.in is missing physical parameter rows.' }
    $impact = $rows[2].Replace(',',' ') -split '\s+'
    $surface = $rows[5].Replace(',',' ') -split '\s+'
    $momentum = $rows[12].Replace(',',' ') -split '\s+'
    $values = @([string]$Cfg.Shards, [string]($Cfg.Events/$Cfg.Shards), [string]$Cfg.LegacyRows[2], [string]$Cfg.LegacyRows[4], [string]$impact[3], [string]$rows[4], (($surface[0..2]) -join ', '), (($surface[3..4]) -join ', '), [string]$momentum[4])
    $labels = @('NFILE','NRUN_PER','NAME','T_PER','X00CUT','alpha,beta,gamma','gsur,gtau,eta','gsym,symm','momp')
    $out = for ($i=0; $i -lt 9; $i++) { $values[$i] + ' ! ' + $labels[$i] }
    [IO.File]::WriteAllLines((Join-Path $Cfg.SourceDirectory 'Input file/INPUT_ANA.in'),[string[]]$out,$script:Utf8NoBom)
}

function Apply-LegacyAnalysisControl {
    param([string]$Path, $Cfg)
    if (-not $Cfg.Legacy -or -not $Cfg.LegacyAnalysisRows) { return }
    $lines = [IO.File]::ReadAllLines($Path)
    Set-FirstValue $lines 'Simulation' 1 $Cfg.LegacyAnalysisRows[4]
    for ($i=0; $i -lt 4; $i++) { Set-FirstValue $lines 'Dynamics' $i $Cfg.LegacyAnalysisRows[5+$i] }
    [IO.File]::WriteAllLines($Path,$lines,$script:Utf8NoBom)
}

function Set-FirstValue {
    param([string[]]$Lines, [string]$Section, [int]$Row, [string]$Values)
    $current = ''
    $found = 0
    for ($i = 0; $i -lt $Lines.Count; $i++) {
        $trim = $Lines[$i].Trim()
        if ($trim -match '^%%\s*(.+)$') { $current = $Matches[1].Trim(); continue }
        if ($current -ne $Section -or -not $trim -or $trim.StartsWith('!') -or $trim.StartsWith('#')) { continue }
        if ($found -eq $Row) {
            $parts = $Lines[$i] -split '!', 2
            $suffix = if ($parts.Count -eq 2) { ' !' + $parts[1] } else { '' }
            $Lines[$i] = $Values + $suffix
            return
        }
        $found++
    }
    throw "Could not update data row $Row in section $Section."
}

function Set-PlainFirstLineValues {
    param([string[]]$Lines, [string]$Regex, [string]$Values)
    for ($i = 0; $i -lt $Lines.Count; $i++) {
        if ($Lines[$i] -match $Regex) {
            $match = [regex]::Match($Lines[$i], $Regex)
            $comment = $Lines[$i].Substring($match.Length)
            $Lines[$i] = $Values + $comment
            return
        }
    }
    throw 'Could not identify required numeric input row.'
}

function Get-Shards {
    param([int]$Events, [int]$Count, [int]$BaseSeed)
    $base = [math]::Floor($Events / $Count)
    $extra = $Events % $Count
    $first = 1
    for ($i = 1; $i -le $Count; $i++) {
        $n = [int]$base + $(if ($i -le $extra) { 1 } else { 0 })
        $last = $first + $n - 1
        $seed = $BaseSeed + $i * 1009
        if ($seed -gt 2000000000) { $seed = 1 + ($seed % 1999999999) }
        [pscustomobject]@{ Index = $i; First = $first; Last = $last; Events = $n; Seed = $seed }
        $first = $last + 1
    }
}

function Get-AttemptDirectory {
    param([string]$Base, [string]$Label)
    $parent = Join-Path $Base $Label
    New-Item -ItemType Directory -Path $parent -Force | Out-Null
    $existing = @(Get-ChildItem -LiteralPath $parent -Directory -Filter 'attempt-*' -ErrorAction SilentlyContinue)
    $n = 1
    if ($existing.Count) {
        $numbers = @($existing | ForEach-Object { if ($_.Name -match '^attempt-(\d+)$') { [int]$Matches[1] } })
        if ($numbers.Count) { $n = ($numbers | Measure-Object -Maximum).Maximum + 1 }
    }
    $path = Join-Path $parent ('attempt-{0:d4}' -f $n)
    New-Item -ItemType Directory -Path $path | Out-Null
    return $path
}

function Set-LqmdInput {
    param([string]$InputPath, $Shard, [int]$Model = 0)
    $lines = [IO.File]::ReadAllLines($InputPath)
    Set-PlainFirstLineValues $lines '^\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+' "$($Shard.First), $($Shard.Last), 0, 0"
    [IO.File]::WriteAllLines($InputPath, $lines, $script:Utf8NoBom)
    $massPath = Join-Path (Split-Path -Parent $InputPath) 'm_select.in'
    $mass = [IO.File]::ReadAllLines($massPath)
    $numericRows = 0
    for ($i = 0; $i -lt $mass.Count; $i++) {
        if ($mass[$i] -match '^\s*([+-]?\d+)\s*(?:!.*)?$') {
            $numericRows++
            $value = switch ($numericRows) { 1 { if ($Model) {$Model} else {0} } 2 { 1 } 3 { $Shard.Seed } default { $null } }
            if ($null -ne $value) { $mass[$i] = [regex]::Replace($mass[$i], '^\s*[+-]?\d+', [string]$value) }
            if ($numericRows -eq 3) { break }
        }
    }
    if ($numericRows -lt 3) { throw 'm_select.in must contain at least three numeric data rows.' }
    [IO.File]::WriteAllLines($massPath, $mass, $script:Utf8NoBom)
}

function Set-AnalysisInput {
    param([string]$InputPath, $Shard)
    $lines = [IO.File]::ReadAllLines($InputPath)
    Set-FirstValue $lines 'CONT_ANA' 0 '0'
    Set-FirstValue $lines 'CONT_ANA' 1 '0'
    Set-FirstValue $lines 'CONT_ANA' 2 "1, $($Shard.Events)"
    Set-FirstValue $lines 'Simulation' 0 ([string]$Shard.Events)
    [IO.File]::WriteAllLines($InputPath, $lines, $script:Utf8NoBom)
    $massPath = Join-Path (Split-Path -Parent $InputPath) 'm_select.in'
    $mass = [IO.File]::ReadAllLines($massPath)
    $numericRows = 0
    for ($i = 0; $i -lt $mass.Count; $i++) {
        if ($mass[$i] -match '^\s*([+-]?\d+)\s*(?:!.*)?$') {
            $numericRows++
            $value = switch ($numericRows) { 1 { 1 } 2 { 1 } 3 { $Shard.Seed } default { $null } }
            if ($null -ne $value) { $mass[$i] = [regex]::Replace($mass[$i], '^\s*[+-]?\d+', [string]$value) }
            if ($numericRows -eq 3) { break }
        }
    }
    if ($numericRows -lt 3) { throw 'm_select.in must contain at least three numeric data rows.' }
    [IO.File]::WriteAllLines($massPath, $mass, $script:Utf8NoBom)
}

function Assert-PhysicsInputs {
    param($Cfg)
    $lqPath = Join-Path $Cfg.SourceDirectory 'Input file/LQMD.in'
    $anaPath = Join-Path $Cfg.SourceDirectory 'Input file/anahyp.in'
    if ($Cfg.Legacy) { $lqPath = Join-Path $Cfg.SourceDirectory 'Input file/LQMD0.in'; $anaPath = Join-Path $Cfg.SourceDirectory 'Input file/anahyp0.in' }
    $lq = Get-Content -LiteralPath $lqPath -Encoding UTF8
    if ($lq[0] -notmatch '^\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+') { throw 'LQMD.in first row is not NRMIN,NRMAX,NTH_RUN,NTH_NOH.' }
    $sel = Get-Content -LiteralPath (Join-Path $Cfg.SourceDirectory 'Input file/m_select.in') -Encoding UTF8
    $nums = @($sel | Where-Object { $_ -match '^\s*\d+\s*(?:!.*)?$' } | ForEach-Object { [int]([regex]::Match($_,'^\s*(\d+)').Groups[1].Value) })
    if ($nums.Count -lt 3 -or $nums[1] -notin @(0,1)) { throw 'm_select.in mass_framework must be 0 or 1; the controller selects framework=1 in its private run copies.' }
    $ana = Get-Content -LiteralPath $anaPath -Encoding UTF8
    $coalLq = $null
    foreach ($line in $lq) { if ($line -match '!\s*SURF_COAL') { $data = ($line -split '!',2)[0].Trim(); if ($data -match '^\s*(\d+)') { $coalLq = [int]$Matches[1]; break } } }
    $section = ''; $coalAna = $null
    foreach ($line in $ana) { if ($line -match '^%%\s*(.+)$') { $section = $Matches[1]; continue }; if ($section -eq 'Models' -and $line -match '!\s*SURF_COAL') { $data = ($line -split '!',2)[0].Trim(); if ($data -match '^\s*(\d+)') { $coalAna = [int]$Matches[1] }; break } }
    if ($null -eq $coalLq -or $null -eq $coalAna) { throw 'Cannot identify both SURF_COAL settings; keep their input comments.' }
    if ($coalLq -ne $coalAna) { throw 'LQMD/anahyp SURF_COAL values differ.' }
    if ($coalLq -ne 0) { throw 'The current controller safely supports SURF_COAL=0 only.' }
    $timeData = Get-SectionData $ana 'Simulation' 2
    if (($timeData -replace '[,\s]+', ',').Trim(',') -ne '1,1,1') { throw 'This controller requires kktot,kkmin,kkmax=1,1,1; it will not change the time-slice physics settings.' }
    $yieldData = Get-SectionData $ana 'Yields' 0
    if (($yieldData -split '[,\s]+')[0] -ne '1') { throw 'CNPCS output must be enabled in the source anahyp.in.' }
    $lqData = (($lq[2] -split '!',2)[0].Trim().Replace(',', ' ')) -split '\s+'
    if ($lqData.Count -lt 4) { throw 'LQMD.in must provide x00cut on row 3.' }
    $lqBcut = [double]::Parse($lqData[3].Replace('D','E').Replace('d','e'), [Globalization.CultureInfo]::InvariantCulture)
    if ($lqBcut -le 0 -or $lqData.Count -lt 5 -or [int]$lqData[4] -ne 1) { throw 'Use a positive explicit x00cut and nx00=1 for one normalized event ensemble.' }
    $anaBcut = Get-Bcut $anaPath
    if ([Math]::Abs($lqBcut - $anaBcut) -gt 0.000001) { throw "LQMD x00cut $lqBcut and anahyp X00CUT $anaBcut differ." }
    $Cfg | Add-Member -NotePropertyName Bcut -NotePropertyValue $lqBcut -Force
    $species = @(Get-WindowsSpecies (Join-Path $Cfg.SourceDirectory 'Input file/CNPCS_in.in'))
    foreach ($target in $script:targets) {
        if (@($species | Where-Object { $_.Z -eq $target.Z -and $_.N -eq $target.N }).Count -ne 1) { throw "CNPCS_in.in must select $($target.Name) exactly once." }
    }
}

function Get-SectionData {
    param([string[]]$Lines, [string]$Section, [int]$Row)
    $current = ''; $index = 0
    foreach ($line in $Lines) {
        if ($line -match '^\s*%%\s*(.+)$') { $current=$Matches[1].Trim(); $index=0; continue }
        $data = ($line -split '!',2)[0].Trim()
        if ($current -ne $Section -or -not $data -or $data.StartsWith('#')) { continue }
        if ($index -eq $Row) { return $data }
        $index++
    }
    throw "Cannot identify $Section data row $Row."
}

function Get-Bcut {
    param([string]$Path)
    $lines = [IO.File]::ReadAllLines($Path)
    $section = ''
    $row = 0
    foreach ($line in $lines) {
        if ($line -match '^%%\s*(.+)$') { $section=$Matches[1].Trim(); $row=0; continue }
        if ($section -ne 'Simulation') { continue }
        $data = ($line -split '!',2)[0].Trim()
        if (-not $data -or $data.StartsWith('#')) { continue }
        if ($row -eq 1) {
            $token = (($data -replace ',',' ').Trim() -split '\s+')[0]
            return [double]::Parse($token.Replace('D','E').Replace('d','e'), [Globalization.CultureInfo]::InvariantCulture)
        }
        $row++
    }
    throw "Could not identify X00CUT in $Path"
}

function New-LqmdTask {
    param($Cfg, $Shard, [switch]$PrepareOnly)
    $label = 'shard-{0:d4}' -f $Shard.Index
    $base = Join-Path $Cfg.RunDirectory 'kinematics'
    if (@(Get-CurrentAttempts $Cfg.RunDirectory 'kinematics' $label).Count) { return $null }
    $prepared = Get-PreparedTask $base $label
    if ($prepared) { return [pscustomobject]@{ Path=$prepared; Task=(Join-Path $prepared 'task.json'); Shard=$Shard } }
    $attempt = Get-AttemptDirectory $base $label
    Copy-RunInputs $Cfg.SourceDirectory $attempt -UseTemplates:$Cfg.Legacy
    Set-LqmdInput (Join-Path $attempt 'Input file/LQMD.in') $Shard $(if ($Cfg.Legacy) {$Cfg.SelectedModel} else {0})
    $state = [ordered]@{ stage='L'; shard=$Shard.Index; first_event=$Shard.First; last_event=$Shard.Last; events=$Shard.Events; seed=$Shard.Seed; status='prepared'; attempt=(Split-Path -Leaf $attempt); created_utc=[DateTime]::UtcNow.ToString('o'); input_hashes=(Get-InputSnapshot $attempt -SkipInitialFiles) }
    Write-JsonAtomic (Join-Path $attempt 'state.json') $state
    $taskPath = Join-Path $attempt 'task.json'
    Write-JsonAtomic $taskPath ([ordered]@{ Executable=$Cfg.LqmdExecutable; ExecutableHash=(Get-Sha256 $Cfg.LqmdExecutable); WorkingDirectory=$attempt; Marker='The program end normally, ok!'; RequiredFiles=@('coll.dat','initangm.dat'); StateFile=(Join-Path $attempt 'state.json') })
    return [pscustomobject]@{ Path=$attempt; Task=$taskPath; Shard=$Shard }
}

function New-AnalysisTask {
    param($Cfg, $Shard, [int]$Model, [string]$HeaderPath, [switch]$PrepareOnly)
    $label = 'shard-{0:d4}' -f $Shard.Index
    $base = Join-Path (Join-Path $Cfg.RunDirectory 'analysis') ('model-{0}' -f $Model)
    if (@(Get-CurrentAttempts $Cfg.RunDirectory ('analysis/model-{0}' -f $Model) $label).Count) { return $null }
    $prepared = Get-PreparedTask $base $label
    if ($prepared) { return [pscustomobject]@{ Path=$prepared; Task=(Join-Path $prepared 'task.json'); Shard=$Shard; Model=$Model } }
    $attempt = Get-AttemptDirectory $base $label
    Copy-RunInputs $Cfg.SourceDirectory $attempt -UseTemplates:$Cfg.Legacy
    $kinAttempts = @(Get-CurrentAttempts $Cfg.RunDirectory 'kinematics' $label)
    if (-not $kinAttempts.Count) { throw "Missing LQMD shard: $label" }
    $kin = $kinAttempts[0].Path
    $collSource = Join-Path $kin 'coll.dat'
    $collTarget = Join-Path $attempt 'coll.dat'
    if ($Shard.Index -eq 1) { Copy-Item -LiteralPath $collSource -Destination $collTarget }
    else {
        $header = [IO.File]::ReadAllLines($HeaderPath)
        if ($header.Count -ne 3 -or $header[0].Trim() -ne '') { throw 'Invalid coll.dat header template; expected three lines including the initial blank record.' }
        $raw = [IO.File]::ReadAllLines($collSource)
        if ($raw.Count -lt 1 -or $raw[0].Trim() -eq '') { throw "Shard $($Shard.Index) unexpectedly contains a coll.dat header." }
        $all = New-Object 'System.Collections.Generic.List[string]'
        $header | ForEach-Object { $all.Add($_) }
        $raw | ForEach-Object { $all.Add($_) }
        [IO.File]::WriteAllLines($collTarget, $all, $script:Utf8NoBom)
    }
    Copy-Item -LiteralPath (Join-Path $kin 'initangm.dat') -Destination (Join-Path $attempt 'initangm.dat')
    $shardA = [pscustomobject]@{ Index=$Shard.Index; First=1; Last=$Shard.Events; Events=$Shard.Events; Seed=$Shard.Seed }
    Set-AnalysisInput (Join-Path $attempt 'Input file/anahyp.in') $shardA
    Apply-LegacyAnalysisControl (Join-Path $attempt 'Input file/anahyp.in') $Cfg
    $selPath = Join-Path $attempt 'Input file/m_select.in'
    $sel = [IO.File]::ReadAllLines($selPath); $numeric = 0
    for ($i=0; $i -lt $sel.Count; $i++) { if ($sel[$i] -match '^\s*([+-]?\d+)\s*(?:!.*)?$') { $numeric++; if ($numeric -eq 1) { $sel[$i] = [regex]::Replace($sel[$i], '^\s*[+-]?\d+', [string]$Model); break } } }
    [IO.File]::WriteAllLines($selPath, $sel, $script:Utf8NoBom)
    $state = [ordered]@{ stage='A'; model=$Model; shard=$Shard.Index; global_first_event=$Shard.First; global_last_event=$Shard.Last; events=$Shard.Events; seed=$Shard.Seed; status='prepared'; attempt=(Split-Path -Leaf $attempt); created_utc=[DateTime]::UtcNow.ToString('o'); input_hashes=(Get-InputSnapshot $attempt -SkipInitialFiles) }
    Write-JsonAtomic (Join-Path $attempt 'state.json') $state
    $taskPath = Join-Path $attempt 'task.json'
    Write-JsonAtomic $taskPath ([ordered]@{ Executable=$Cfg.ExfragExecutable; ExecutableHash=(Get-Sha256 $Cfg.ExfragExecutable); WorkingDirectory=$attempt; Marker='OK, the program ends normally !!'; RequiredFiles=@('CNPCS_HOT.dat','CNPCS_COLD.dat'); StateFile=(Join-Path $attempt 'state.json') })
    return [pscustomobject]@{ Path=$attempt; Task=$taskPath; Shard=$Shard; Model=$Model }
}

function Invoke-Tasks {
    param([object[]]$Tasks, [int]$Limit, [string]$ConfigPath)
    $queue = New-Object 'System.Collections.Generic.Queue[object]'
    foreach ($task in $Tasks) { $queue.Enqueue($task) }
    $running = New-Object 'System.Collections.Generic.List[object]'
    $powershell = (Get-Command $script:ChildPowerShell -ErrorAction Stop).Source
    while ($queue.Count -gt 0 -or $running.Count -gt 0) {
        while ($queue.Count -gt 0 -and $running.Count -lt $Limit) {
            $task = $queue.Dequeue()
            Write-Output "Starting $Stage task: $($task.Task)"
            $out = Join-Path $task.Path 'controller.stdout.log'
            $err = Join-Path $task.Path 'controller.stderr.log'
            $args = @('-NoProfile','-ExecutionPolicy','Bypass','-File', (Quote-ProcessArgument $PSCommandPath),'-Stage', $Stage,'-Config', (Quote-ProcessArgument $ConfigPath),'-InternalTask', (Quote-ProcessArgument $task.Task)) -join ' '
            $proc = Start-ControllerProcess @{FilePath=$powershell;ArgumentList=$args;WorkingDirectory=$PSScriptRoot;RedirectStandardOutput=$out;RedirectStandardError=$err}
            $null = $proc.Handle
            $running.Add([pscustomobject]@{ Process=$proc; Task=$task })
        }
        Start-Sleep -Milliseconds 500
        foreach ($entry in $running.ToArray()) {
            $entry.Process.Refresh()
            if ($entry.Process.HasExited) {
                if ($entry.Process.ExitCode -ne 0) { throw "Worker failed (exit $($entry.Process.ExitCode)): $($entry.Task.Path)" }
                $statePath = Join-Path $entry.Task.Path 'state.json'
                $state = Read-JsonShared $statePath
                if ($state.status -ne 'completed') { throw "Worker ended without a completed state: $($entry.Task.Path)" }
                Write-Output "Completed $Stage task: $($entry.Task.Task)"
                $running.Remove($entry) | Out-Null
                $entry.Process.Dispose()
            }
        }
    }
}

function Quote-ProcessArgument {
    param([string]$Value)
    if ($Value.Contains('"')) { throw 'Process arguments may not contain a quote character.' }
    return '"' + $Value + '"'
}

function Start-ControllerProcess {
    param([hashtable]$Parameters)
    $Parameters.PassThru = $true
    if ($script:WindowsHost) { $Parameters.WindowStyle = 'Hidden' }
    Start-Process @Parameters
}

function Get-ComputeSlot {
    param([switch]$NoWait)
    if ($script:WindowsHost) {
        $slot = New-Object System.Threading.Semaphore(25,25,'Local\SpallationWindowsComputeSlots25')
        if ($NoWait) { if (-not $slot.WaitOne(0)) { $slot.Dispose(); return $null } }
        else { $slot.WaitOne() | Out-Null }
        return $slot
    }
    $directory = Join-Path ([Environment]::GetFolderPath('UserProfile')) '.local/state/spallation/compute-slots'
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
    do {
        foreach ($index in 1..25) {
            if ($script:HeldSlots.ContainsKey($index)) { continue }
            $slot = [IO.File]::Open((Join-Path $directory "$index.lock"),[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::ReadWrite)
            try { $slot.Lock(0,1); $script:HeldSlots[$index]=$true; return $slot }
            catch [IO.IOException] { $slot.Dispose() }
        }
        if ($NoWait) { return $null }
        Start-Sleep -Milliseconds 250
    } while ($true)
}

function Close-ComputeSlot {
    param($Slot)
    if (-not $Slot) { return }
    try {
        if ($script:WindowsHost) { $Slot.Release() | Out-Null }
        else {
            $Slot.Unlock(0,1)
            $script:HeldSlots.Remove([int][IO.Path]::GetFileNameWithoutExtension($Slot.Name))
        }
    } finally { $Slot.Dispose() }
}

function Enter-ControllerMutex {
    $name = 'Local\SpallationWindowsRunController'
    if (-not $script:WindowsHost) { $name = 'SpallationRunController' }
    if ($script:ConfigIsLegacy -and $script:ConfigRunDirectory) {
        $identity = $script:ConfigRunDirectory
        if ($script:WindowsHost) { $identity = $identity.ToLowerInvariant() }
        $bytes = [Text.Encoding]::UTF8.GetBytes($identity)
        $sha = [Security.Cryptography.SHA256]::Create()
        try { $suffix = ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-','').Substring(0,24) } finally { $sha.Dispose() }
        $name += '_' + $suffix
    }
    $mutex = New-Object System.Threading.Mutex($false, $name)
    try { if (-not $mutex.WaitOne(0)) { throw 'Another Windows batch controller is active. Wait for it to finish before starting the other stage.' } }
    catch [System.Threading.AbandonedMutexException] { }
    return $mutex
}

function Assert-NoActiveWorkers {
    param($Cfg)
    if ($Cfg.Legacy) { return }
    $active = @(Get-Process -Name LQMD,Console1,exfrag -ErrorAction SilentlyContinue)
    if ($active.Count -gt 0) { throw 'Existing LQMD/exfrag/Console1 processes are running. Wait for them to exit before starting a batch.' }
    $taken = New-Object 'System.Collections.Generic.List[object]'
    try {
        while ($taken.Count -lt 25) {
            $slot = Get-ComputeSlot -NoWait
            if (-not $slot) { break }
            $taken.Add($slot)
        }
        if ($taken.Count -ne 25) { throw 'Prior controller workers are still active. Wait for their shared compute slots to be released.' }
    } finally {
        foreach ($slot in $taken) { Close-ComputeSlot $slot }
    }
}

function Assert-LegacyRunIdle {
    param([string]$RunDirectory)
    if (-not (Test-Path -LiteralPath $RunDirectory -PathType Container)) { return }
    foreach ($file in Get-ChildItem -LiteralPath $RunDirectory -Filter state.json -File -Recurse -ErrorAction SilentlyContinue) {
        $state = Read-JsonShared $file.FullName
        if ($state.status -ne 'running') { continue }
        foreach ($id in @($state.program_pid, $state.pid)) {
            if ($id) {
                $process = Get-Process -Id ([int]$id) -ErrorAction SilentlyContinue
                if ($process) { throw "A worker for this run is still active (PID $id, state $($file.FullName))." }
            }
        }
    }
}

function Invoke-Worker {
    param([string]$TaskPath)
    $task = Read-JsonShared $TaskPath
    $statePath = $task.StateFile
    $state = Read-JsonShared $statePath
    if ($state.status -ne 'prepared') { throw "Task state is not prepared: $statePath" }
    Assert-PrivateInputs $task.WorkingDirectory $state
    if ((Get-Sha256 $task.Executable) -cne $task.ExecutableHash) { throw 'Prepared executable has changed.' }
    $state.status = 'running'
    $state | Add-Member -NotePropertyName started_utc -NotePropertyValue ([DateTime]::UtcNow.ToString('o')) -Force
    $state | Add-Member -NotePropertyName pid -NotePropertyValue $PID -Force
    Set-JsonAtomic $statePath $state
    $slot = $null
    try {
        $slot = Get-ComputeSlot
        foreach ($envName in @('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS')) { [Environment]::SetEnvironmentVariable($envName, '1', 'Process') }
        $proc = Start-ControllerProcess @{FilePath=$task.Executable;WorkingDirectory=$task.WorkingDirectory;RedirectStandardOutput=(Join-Path $task.WorkingDirectory 'program.stdout.log');RedirectStandardError=(Join-Path $task.WorkingDirectory 'program.stderr.log')}
        $null = $proc.Handle
        $state | Add-Member -NotePropertyName program_pid -NotePropertyValue $proc.Id -Force
        Set-JsonAtomic $statePath $state
        try { $proc.PriorityClass = [Diagnostics.ProcessPriorityClass]::BelowNormal } catch { }
        $proc.WaitForExit()
        $proc.Refresh()
        if ($proc.ExitCode -ne 0) { throw "Program exit code $($proc.ExitCode)." }
        $stdout = Join-Path $task.WorkingDirectory 'program.stdout.log'
        $stderr = Join-Path $task.WorkingDirectory 'program.stderr.log'
        $text = ''
        if (Test-Path -LiteralPath $stdout) { $text += [IO.File]::ReadAllText($stdout) }
        if (Test-Path -LiteralPath $stderr) { $text += [IO.File]::ReadAllText($stderr) }
        if (-not $text.Contains([string]$task.Marker)) { throw "Missing normal-completion marker: $($task.Marker)" }
        if ($text -match '(?im)FAULT\s*!|fortran runtime error|end of file|segmentation fault|fatal error') { throw 'Program log contains a Fortran/runtime fatal error marker.' }
        foreach ($file in $task.RequiredFiles) {
            $path = Join-Path $task.WorkingDirectory $file
            if (-not (Test-Path -LiteralPath $path -PathType Leaf) -or (Get-Item -LiteralPath $path).Length -eq 0) { throw "Missing or empty output: $file" }
        }
        if ($state.stage -eq 'L') {
            $initLines = [IO.File]::ReadAllLines((Join-Path $task.WorkingDirectory 'initangm.dat')).Count
            if ($initLines -ne [int]$state.events) { throw "Expected $($state.events) initangm records; found $initLines." }
        } else {
            foreach ($name in @('CNPCS_HOT.dat','CNPCS_COLD.dat')) { $null = Read-WindowsCnpcs (Join-Path $task.WorkingDirectory $name) }
        }
        $outputHashes = [ordered]@{}
        Assert-PrivateInputs $task.WorkingDirectory $state
        foreach ($file in $task.RequiredFiles) { $outputHashes[$file] = Get-Sha256 (Join-Path $task.WorkingDirectory $file) }
        $state.status = 'completed'
        $state | Add-Member -NotePropertyName exit_code -NotePropertyValue ([int]$proc.ExitCode) -Force
        $state | Add-Member -NotePropertyName output_hashes -NotePropertyValue $outputHashes -Force
        $state | Add-Member -NotePropertyName finished_utc -NotePropertyValue ([DateTime]::UtcNow.ToString('o')) -Force
        Set-JsonAtomic $statePath $state
    } catch {
        $state.status = 'failed'
        $state | Add-Member -NotePropertyName error -NotePropertyValue $_.Exception.Message -Force
        $state | Add-Member -NotePropertyName finished_utc -NotePropertyValue ([DateTime]::UtcNow.ToString('o')) -Force
        Set-JsonAtomic $statePath $state
        throw
    } finally {
        Close-ComputeSlot $slot
    }
}

function Get-StateFiles {
    param([string]$Root, [string]$StageName, [int]$Model = 0)
    $base = Join-Path $Root 'kinematics'
    if ($StageName -eq 'A') { $base = Join-Path (Join-Path $Root 'analysis') ('model-{0}' -f $Model) }
    $states = @()
    Get-ChildItem -LiteralPath $base -Filter state.json -File -Recurse -ErrorAction SilentlyContinue | ForEach-Object {
        $state = Get-Content -LiteralPath $_.FullName -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($state.stage -eq $StageName -and ($StageName -ne 'A' -or [int]$state.model -eq $Model)) { $states += $state }
    }
    return ,$states
}

function Get-CurrentAttempts {
    param([string]$Root, [string]$Parent, [string]$Label)
    $dir = Join-Path (Join-Path $Root $Parent) $Label
    if (-not (Test-Path -LiteralPath $dir)) { return @() }
    $items = @(Get-ChildItem -LiteralPath $dir -Directory -Filter 'attempt-*' | Sort-Object Name -Descending)
    if (-not $items.Count) { return @() }
    foreach ($item in $items) {
    $stateFile = Join-Path $item.FullName 'state.json'
    if (-not (Test-Path -LiteralPath $stateFile)) { continue }
    $state = Read-JsonShared $stateFile
    if ($state.status -eq 'completed') {
        Assert-PrivateInputs $item.FullName $state
        if (-not $state.output_hashes) { throw "Completed shard has no output hashes: $stateFile" }
        foreach ($name in $state.output_hashes.PSObject.Properties.Name) {
            $output = Join-Path $item.FullName $name
            if (-not (Test-Path -LiteralPath $output -PathType Leaf)) { throw "Completed output is missing: $output" }
            if ((Get-Sha256 $output) -cne [string]$state.output_hashes.$name) { throw "Completed output was modified: $output" }
        }
        return @([pscustomobject]@{ Path=$item.FullName; State=$state; Done=$true })
    }
    }
    return @()
}

function Assert-PrivateInputs {
    param([string]$Directory, $State)
    if (-not $State.input_hashes) { throw "Missing private input snapshot: $Directory" }
    $now = Get-InputSnapshot $Directory -SkipInitialFiles
    if (($now | ConvertTo-Json -Compress -Depth 10) -cne ($State.input_hashes | ConvertTo-Json -Compress -Depth 10)) { throw "Private run inputs changed: $Directory" }
}

function Get-PreparedTask {
    param([string]$Base, [string]$Label)
    $parent = Join-Path $Base $Label
    if (-not (Test-Path -LiteralPath $parent)) { return $null }
    $items = @(Get-ChildItem -LiteralPath $parent -Directory -Filter 'attempt-*' | Sort-Object Name -Descending)
    if (-not $items.Count) { return $null }
    $stateFile = Join-Path $items[0].FullName 'state.json'
    $taskFile = Join-Path $items[0].FullName 'task.json'
    if (-not (Test-Path -LiteralPath $stateFile) -or -not (Test-Path -LiteralPath $taskFile)) { return $null }
    $state = Read-JsonShared $stateFile
    if ($state.status -eq 'prepared') { Assert-PrivateInputs $items[0].FullName $state; return $items[0].FullName }
    return $null
}

function Merge-Lqmd {
    param($Cfg, [object[]]$Shards)
    $dest = Join-Path (Join-Path $Cfg.RunDirectory 'kinematics') 'merged'
    if (Test-Path -LiteralPath $dest) {
        foreach ($name in @('coll.dat','initangm.dat','merge.json')) { if (-not (Test-Path -LiteralPath (Join-Path $dest $name) -PathType Leaf)) { throw "Incomplete existing LQMD merge: $dest" } }
        return
    }
    $temp = $dest + '.tmp-' + [guid]::NewGuid().ToString('N')
    New-Item -ItemType Directory -Path $temp | Out-Null
    foreach ($name in @('coll.dat','initangm.dat')) {
        $outPath = Join-Path $temp $name
        $out = [IO.File]::Open($outPath, [IO.FileMode]::CreateNew, [IO.FileAccess]::Write, [IO.FileShare]::None)
        try {
            foreach ($shard in $Shards) {
                $attempt = @(Get-CurrentAttempts $Cfg.RunDirectory 'kinematics' ('shard-{0:d4}' -f $shard.Index))
                if (-not $attempt.Count) { throw "LQMD shard $($shard.Index) is incomplete or changed." }
                $source = [IO.File]::OpenRead((Join-Path $attempt[0].Path $name))
                try { $source.CopyTo($out) } finally { $source.Dispose() }
            }
        } finally { $out.Dispose() }
    }
    $meta = [ordered]@{ events=$Cfg.Events; shards=$Cfg.Shards; coll_sha256=(Get-Sha256 (Join-Path $temp 'coll.dat')); initangm_sha256=(Get-Sha256 (Join-Path $temp 'initangm.dat')); created_utc=[DateTime]::UtcNow.ToString('o') }
    Write-JsonAtomic (Join-Path $temp 'merge.json') $meta
    [IO.Directory]::Move($temp,$dest)
}

function Assert-AllComplete {
    param($Cfg, [object[]]$Shards, [string]$StageName, [int]$Model = 0)
    foreach ($shard in $Shards) {
        $label = 'shard-{0:d4}' -f $shard.Index
        $parent = if ($StageName -eq 'L') { 'kinematics' } else { 'analysis/model-{0}' -f $Model }
        $done = @(Get-CurrentAttempts $Cfg.RunDirectory $parent $label)
        if (-not $done.Count -or $done[0].State.status -ne 'completed') { throw "$StageName is missing a valid completed shard: $parent/$label" }
        if ([int]$done[0].State.events -ne [int]$shard.Events) { throw "Shard event count does not match plan: $parent/$label" }
    }
}

function Set-MergedModelInputs {
    param([string]$InputDirectory, $Cfg, [int]$Model)
    $anaPath = Join-Path $InputDirectory 'anahyp.in'
    $shard = [pscustomobject]@{ Events=$Cfg.Events; Seed=$Cfg.BaseSeed }
    Set-AnalysisInput $anaPath $shard
    $massPath = Join-Path $InputDirectory 'm_select.in'
    $mass = [IO.File]::ReadAllLines($massPath)
    $numeric = 0
    for ($i=0; $i -lt $mass.Count; $i++) {
        if ($mass[$i] -match '^\s*[+-]?\d+\s*(?:!.*)?$') {
            $numeric++
            if ($numeric -eq 1) { $mass[$i] = [regex]::Replace($mass[$i], '^\s*[+-]?\d+', [string]$Model); break }
        }
    }
    [IO.File]::WriteAllLines($massPath, $mass, $script:Utf8NoBom)
}

function Merge-Cnpcs {
    param($Cfg, [object[]]$Shards, [int]$Model)
    $dest = Join-Path (Join-Path $Cfg.RunDirectory 'merged') ('model-{0}' -f $Model)
    if (Test-Path -LiteralPath $dest) {
        $oldMetaPath = Join-Path $dest 'merge.json'
        if (-not (Test-Path -LiteralPath $oldMetaPath)) { throw "Refusing to overwrite incomplete merged output: $dest" }
        $oldMeta = Get-Content -LiteralPath $oldMetaPath -Raw -Encoding UTF8 | ConvertFrom-Json
        if ([int]$oldMeta.events -ne $Cfg.Events -or [int]$oldMeta.model -ne $Model) { throw "Merged output metadata mismatch: $dest" }
        foreach ($name in @('CNPCS_HOT.dat','CNPCS_COLD.dat')) {
            if (-not (Test-Path -LiteralPath (Join-Path $dest $name)) -or (Get-Sha256 (Join-Path $dest $name)) -cne [string]$oldMeta.output_hashes.$name) { throw "Existing merged output was modified: $dest/$name" }
        }
        return
    }
    $temp = $dest + '.tmp-' + [guid]::NewGuid().ToString('N')
    New-Item -ItemType Directory -Path $temp -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $Cfg.SourceDirectory 'Input file') -Destination $temp -Recurse
    Set-MergedModelInputs (Join-Path $temp 'Input file') $Cfg $Model
    foreach ($name in @('CNPCS_HOT.dat','CNPCS_COLD.dat')) {
        $sums = @{ Z=@{}; N=@{}; A=@{}; I=@{} }
        $template = $null
        $totalEvents = 0
        foreach ($shard in $Shards) {
            $label = 'shard-{0:d4}' -f $shard.Index
            $attempt = @(Get-CurrentAttempts $Cfg.RunDirectory ('analysis/model-{0}' -f $Model) $label)
            if (-not $attempt.Count) { throw "No complete analysis attempt for model $Model shard $($shard.Index)." }
            $path = Join-Path $attempt[0].Path $name
            $lines = [IO.File]::ReadAllLines($path)
            if (-not $template) { $template = $lines.Clone() }
            $section = ''
            foreach ($line in $lines) {
                $trim = $line.Trim()
                if ($trim -match 'NUCLIDE\s+CROSS\s+SECTION\s+AS\s+A\s+FUNCTION\s+OF\s+([ZNA])') { $section=$Matches[1]; continue }
                if ($trim -match '^CROSS SECTION OF ISOTOPICS') { $section='I'; continue }
                if (-not $trim) { $section=''; continue }
                if ($section -and $trim -match '^\s*(\d+)\s+([+-]?[0-9.EeDd+-]+)') {
                    $id = [string]$Matches[1]
                    $value = [double]::Parse($Matches[2].Replace('D','E').Replace('d','e'), [Globalization.CultureInfo]::InvariantCulture)
                    if (-not $sums[$section].ContainsKey($id)) { $sums[$section][$id] = [double]0 }
                    $sums[$section][$id] += $value * [int]$shard.Events
                }
            }
            $totalEvents += [int]$shard.Events
        }
        if ($totalEvents -ne $Cfg.Events) { throw "Merged event total mismatch: $totalEvents != $($Cfg.Events)." }
        $speciesCount = 0
        if ($name -eq 'CNPCS_COLD.dat') {
            $speciesLines = [IO.File]::ReadAllLines((Join-Path $temp 'Input file/CNPCS_in.in'))
            for ($i=2; $i -lt $speciesLines.Count; $i++) { $line=$speciesLines[$i].Trim(); if (-not $line -or $line -match '^-{3,}') { break }; $speciesCount++ }
            if ($speciesCount -lt 1) { throw 'CNPCS_in.in has no selected nuclides.' }
        }
        $output = New-Object 'System.Collections.Generic.List[string]'
        for ($i=0; $i -lt $template.Count; $i++) {
            $trim = $template[$i].Trim()
            $section = ''
            if ($trim -match 'NUCLIDE\s+CROSS\s+SECTION\s+AS\s+A\s+FUNCTION\s+OF\s+([ZNA])') { $section=$Matches[1] }
            elseif ($trim -match '^CROSS SECTION OF ISOTOPICS') { $section='I' }
            if ($section) {
                $output.Add($template[$i])
                $max = switch ($section) { 'Z' { 150 } 'N' { 150 } 'A' { 300 } 'I' { $speciesCount } }
                for ($id=1; $id -le $max; $id++) {
                    $value = 0.0
                    if ($sums[$section].ContainsKey([string]$id)) { $value = $sums[$section][[string]$id] / $totalEvents }
                    $output.Add(('{0,6} {1}' -f $id, $value.ToString('E12', [Globalization.CultureInfo]::InvariantCulture)))
                }
                $output.Add('')
                while ($i + 1 -lt $template.Count -and $template[$i+1].Trim()) { $i++ }
            } elseif (-not $trim -or $trim -match '^[\s-]+$') { $output.Add($template[$i]) }
        }
        [IO.File]::WriteAllLines((Join-Path $temp $name), $output, $script:Utf8NoBom)
    }
    $meta = [ordered]@{ model=$Model; model_name=$script:Models[[string]$Model]; events=$Cfg.Events; bcut_fm=$Cfg.Bcut; mass_framework=1; normalization='sum(slice cross section * slice events) / total events'; input_precision_note='Each slice starts from the existing rounded CNPCS text output.'; merged_input_notice='The adjusted Input file is for result provenance/data parsing only; the merged directory is not a runnable simulation restart.'; created_utc=[DateTime]::UtcNow.ToString('o') }
    $meta.output_hashes = [ordered]@{
        'CNPCS_HOT.dat' = Get-Sha256 (Join-Path $temp 'CNPCS_HOT.dat')
        'CNPCS_COLD.dat' = Get-Sha256 (Join-Path $temp 'CNPCS_COLD.dat')
    }
    Write-JsonAtomic (Join-Path $temp 'merge.json') $meta
    [IO.Directory]::Move($temp, $dest)
}

function Assert-RunManifest {
    param($Cfg, [string]$ManifestPath, [object]$Snapshot)
    if (Test-Path -LiteralPath $ManifestPath) {
        $old = Read-JsonShared $ManifestPath
        $expected = [ordered]@{ source=$Cfg.SourceDirectory; events=$Cfg.Events; shards=$Cfg.Shards; max_workers=$Cfg.MaxWorkers; base_seed=$Cfg.BaseSeed; models=@($Cfg.Models); input_hashes=$Snapshot }
        $oldCheck = [ordered]@{ source=$old.source; events=$old.events; shards=$old.shards; max_workers=$old.max_workers; base_seed=$old.base_seed; models=@($old.models); input_hashes=$old.input_hashes }
        if (($expected | ConvertTo-Json -Depth 10 -Compress) -cne ($oldCheck | ConvertTo-Json -Depth 10 -Compress)) { throw 'Run configuration or source hashes changed; use a new RunDirectory.' }
        return $old
    }
    return $null
}

if (-not $Stage) { throw 'Use -Stage L or A.' }
if ($InternalTask) {
    if (-not $Stage) { throw 'InternalTask requires Stage.' }
    Invoke-Worker $InternalTask
    exit 0
}
if (-not $Config) { throw 'Use -Config with INPUT.in, INPUT_ANA.in, or a legacy JSON configuration path.' }
if ($PrepareOnly -and $MergeOnly) { throw '-PrepareOnly and -MergeOnly cannot be combined.' }
$script:Stage = $Stage
$cfg = Get-Config $Config
$autoAnalysis = $null
$script:ConfigIsLegacy = [bool]$cfg.Legacy
$script:ConfigRunDirectory = [string]$cfg.RunDirectory
$mutex = Enter-ControllerMutex
try {
Assert-NoActiveWorkers $cfg
$run = $cfg.RunDirectory
if ($cfg.Legacy) { Assert-LegacyRunIdle $run }
Write-Output "Stage $Stage; events=$($cfg.Events); shards=$($cfg.Shards); workers=$($cfg.MaxWorkers); run=$run"
if (-not (Test-Path -LiteralPath $run)) { New-Item -ItemType Directory -Path $run -Force | Out-Null }
Assert-PhysicsInputs $cfg
$hashes = Get-InputSnapshot $cfg.SourceDirectory @($cfg.LqmdExecutable,$cfg.ExfragExecutable)
$manifestPath = Join-Path $run 'manifest.json'
$oldManifest = Assert-RunManifest $cfg $manifestPath $hashes
if (-not $oldManifest) {
    if (@(Get-ChildItem -LiteralPath $run -Force).Count -gt 0) { throw "RunDirectory is not empty and has no controller manifest: $run" }
    $shards0 = @(Get-Shards $cfg.Events $cfg.Shards $cfg.BaseSeed)
    $manifest = [ordered]@{ format_version=1; source=$cfg.SourceDirectory; run_directory=$run; events=$cfg.Events; shards=$cfg.Shards; max_workers=$cfg.MaxWorkers; base_seed=$cfg.BaseSeed; models=@($cfg.Models); input_hashes=$hashes; shard_plan=@($shards0); created_utc=[DateTime]::UtcNow.ToString('o'); stages=[ordered]@{ L='pending'; A='pending'; merge='pending' } }
    Write-JsonAtomic $manifestPath $manifest
}
$shards = @(Get-Shards $cfg.Events $cfg.Shards $cfg.BaseSeed)
if ($PrepareOnly) {
    if ($Stage -eq 'L') { foreach ($shard in $shards) { New-LqmdTask $cfg $shard -PrepareOnly | Out-Null } }
    else {
        Assert-AllComplete $cfg $shards 'L'
        $first = @(Get-CurrentAttempts $run 'kinematics' 'shard-0001')
        if (-not $first.Count) { throw 'Complete the L stage first; A preparation needs the validated coll header.' }
        $headPath = Join-Path $first[0].Path 'coll.dat'
        $headerPath = Join-Path $run 'coll-header.tmp'
        [IO.File]::WriteAllLines($headerPath, [IO.File]::ReadAllLines($headPath)[0..2], $script:Utf8NoBom)
        foreach ($model in $cfg.Models) { foreach ($shard in $shards) { New-AnalysisTask $cfg $shard $model $headerPath -PrepareOnly | Out-Null } }
        Remove-Item -LiteralPath $headerPath
    }
    Write-Output "Prepared stage $Stage under $run"
    return
}
    if ($MergeOnly -and $Stage -ne 'A') { throw '-MergeOnly is available for Stage A only.' }
    if ($Stage -eq 'L') {
        $tasks = New-Object 'System.Collections.Generic.List[object]'
        foreach ($shard in $shards) {
            $old = @(Get-CurrentAttempts $run 'kinematics' ('shard-{0:d4}' -f $shard.Index))
            if ($old.Count) { continue }
            $task = New-LqmdTask $cfg $shard
            if ($task) { $tasks.Add($task) }
        }
        Invoke-Tasks -Tasks $tasks.ToArray() -Limit $cfg.MaxWorkers -ConfigPath $cfg.ConfigPath
        Assert-AllComplete $cfg $shards 'L'
        Merge-Lqmd $cfg $shards
        $runState = Read-JsonShared $manifestPath
        $runState.stages.L = 'completed'
        Set-JsonAtomic $manifestPath $runState
        if ($cfg.Legacy -and $cfg.CallAnalysis) { $autoAnalysis=Join-Path $cfg.SourceDirectory 'Input file/INPUT_ANA.in' }
    } elseif ($Stage -eq 'A') {
        Assert-AllComplete $cfg $shards 'L'
        if (-not $MergeOnly) {
            $firstAttempt = @(Get-CurrentAttempts $run 'kinematics' 'shard-0001')
            $headPath = Join-Path $firstAttempt[0].Path 'coll.dat'
            $header = [IO.File]::ReadAllLines($headPath)
            if ($header.Count -lt 3 -or $header[0].Trim() -ne '') { throw 'First LQMD coll.dat does not contain the expected three-line header.' }
            $headerPath = Join-Path $run 'coll-header.tmp'
            [IO.File]::WriteAllLines($headerPath, $header[0..2], $script:Utf8NoBom)
            $tasks = New-Object 'System.Collections.Generic.List[object]'
            foreach ($model in $cfg.Models) { foreach ($shard in $shards) {
                $old = @(Get-CurrentAttempts $run ('analysis/model-{0}' -f $model) ('shard-{0:d4}' -f $shard.Index))
                if ($old.Count) { continue }
                $task = New-AnalysisTask $cfg $shard $model $headerPath
                if ($task) { $tasks.Add($task) }
            } }
            Remove-Item -LiteralPath $headerPath
            Invoke-Tasks -Tasks $tasks.ToArray() -Limit $cfg.MaxWorkers -ConfigPath $cfg.ConfigPath
        }
        foreach ($model in $cfg.Models) { Assert-AllComplete $cfg $shards 'A' $model }
        foreach ($model in $cfg.Models) { Merge-Cnpcs $cfg $shards $model }
        Export-WindowsControllerResults $cfg
        $runState = Read-JsonShared $manifestPath
        $runState.stages.A = 'completed'
        $runState.stages.merge = 'completed'
        Set-JsonAtomic $manifestPath $runState
    }
Write-Output "Completed stage $Stage under $run"
} finally {
    if ($mutex) { $mutex.ReleaseMutex(); $mutex.Dispose() }
}
if ($autoAnalysis) {
    Write-Output 'CALL_ANA=1: starting the matching analysis controller.'
    & $PSCommandPath -Stage A -Config $autoAnalysis
}
