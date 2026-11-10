[CmdletBinding()]
param(
    [string]$Compiler = 'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin\ifx.exe',
    [string]$SetVars = 'C:\Program Files (x86)\Intel\oneAPI\setvars.bat',
    [string]$LqmdMockExecutable,
    [string]$ExfragMockExecutable,
    [switch]$SkipBuild,
    [switch]$SkipConcurrency
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$backend = Join-Path $repo 'scripts/windows_event_controller.ps1'
$fixture = Join-Path $PSScriptRoot 'fixtures/windows_controller_mock.f90'
if (-not (Test-Path -LiteralPath $backend -PathType Leaf)) { throw "Controller backend not found: $backend" }
if (-not (Test-Path -LiteralPath $fixture -PathType Leaf)) { throw "Mock fixture not found: $fixture" }

function Assert-True {
    param([bool]$Condition, [string]$Message)
    if (-not $Condition) { throw "ASSERT FAILED: $Message" }
}

function Write-Utf8 {
    param([string]$Path, [string]$Text)
    [IO.File]::WriteAllText($Path, $Text, [Text.UTF8Encoding]::new($false))
}

function Quote-Argument {
    param([string]$Value)
    return '"' + $Value.Replace('"','\"') + '"'
}

function Invoke-Controller {
    param([string[]]$Arguments)
    $powershell = (Get-Command powershell.exe -ErrorAction Stop).Source
    $stdout = Join-Path $script:testRoot ('invoke-' + [guid]::NewGuid().ToString('N') + '.stdout.log')
    $stderr = Join-Path $script:testRoot ('invoke-' + [guid]::NewGuid().ToString('N') + '.stderr.log')
    $quoted = @($Arguments | ForEach-Object { Quote-Argument $_ })
    $process = Start-Process -FilePath $powershell -ArgumentList (@('-NoProfile','-ExecutionPolicy','Bypass') + $quoted) -PassThru -Wait -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr
    $text = ''
    if (Test-Path -LiteralPath $stdout) { $text += [IO.File]::ReadAllText($stdout) }
    if (Test-Path -LiteralPath $stderr) { $text += [IO.File]::ReadAllText($stderr) }
    return [pscustomobject]@{ ExitCode=$process.ExitCode; Text=$text }
}

function Assert-ControllerRejected {
    param([string[]]$Arguments, [string]$Reason)
    $result = Invoke-Controller $Arguments
    Assert-True ($result.ExitCode -ne 0) "controller should reject $Reason"
    return $result.Text
}

function Invoke-Ifx {
    param([string]$CompilerPath, [string]$SetVarsPath, [string]$ExePath, [string]$SourcePath)
    foreach ($item in @($CompilerPath,$SetVarsPath,$SourcePath)) {
        if ($item.Contains('"')) { throw "Compiler path contains an unsupported quote: $item" }
    }
    $command = 'call "' + $SetVarsPath + '" >nul && "' + $CompilerPath + '" /nologo /exe:"' + $ExePath + '" "' + $SourcePath + '"'
    & $env:ComSpec /d /s /c $command
    if ($LASTEXITCODE -ne 0) { throw "ifx fixture build failed with exit code $LASTEXITCODE." }
}

function New-RunConfig {
    param([string]$Path, [string]$RunDirectory, [int]$Events=7, [int]$Shards=3, [int]$Workers=2)
    $value = [ordered]@{
        SourceDirectory = (Join-Path $script:testRoot 'source')
        RunDirectory = $RunDirectory
        LqmdExecutable = $script:lqmdExe
        ExfragExecutable = $script:exfragExe
        Events = $Events
        Shards = $Shards
        MaxWorkers = $Workers
        BaseSeed = 314159
        Models = @(1,2,3,4,5)
    }
    Write-Utf8 $Path ($value | ConvertTo-Json -Depth 4)
    return $value
}

function Read-CnpcsValue {
    param([string]$Path, [string]$Section, [int]$Key)
    $lines = [IO.File]::ReadAllLines($Path)
    $active = $false
    foreach ($line in $lines) {
        $trim = $line.Trim()
        if ($trim -match 'NUCLIDE\s+CROSS\s+SECTION\s+AS\s+A\s+FUNCTION\s+OF\s+([ZNA])') {
            $active = $Matches[1] -eq $Section
            continue
        }
        if ($trim -match '^CROSS SECTION OF ISOTOPICS') { $active = $Section -eq 'I'; continue }
        if (-not $trim) { $active = $false; continue }
        if ($active -and $trim -match '^\s*(\d+)\s+([+-]?[0-9.EeDd+-]+)') {
            if ([int]$Matches[1] -eq $Key) {
                return [double]::Parse($Matches[2].Replace('D','E').Replace('d','e'), [Globalization.CultureInfo]::InvariantCulture)
            }
        }
    }
    throw "Could not find $Section key $Key in $Path"
}

function Assert-Near {
    param([double]$Actual, [double]$Expected, [double]$Tolerance, [string]$Message)
    if ([math]::Abs($Actual - $Expected) -gt $Tolerance) { throw "ASSERT FAILED: $Message (actual=$Actual expected=$Expected)" }
}

function Clear-MergedOutputs {
    param([string]$RunDirectory)
    $rootPath = [IO.Path]::GetFullPath($script:testRoot).TrimEnd('\') + '\'
    $targetPath = [IO.Path]::GetFullPath((Join-Path $RunDirectory 'merged'))
    if (-not $targetPath.StartsWith($rootPath, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clear a merged directory outside the disposable test root: $targetPath"
    }
    if (Test-Path -LiteralPath $targetPath) { Remove-Item -LiteralPath $targetPath -Recurse -Force }
}

$script:testRoot = Join-Path $repo ('tmp/tests/windows-controller/' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $script:testRoot -Force | Out-Null

if (-not $SkipBuild) {
    if (-not (Test-Path -LiteralPath $Compiler -PathType Leaf)) { throw "Intel ifx compiler not found: $Compiler; pass -Compiler or provide prebuilt mocks with -SkipBuild." }
    if (-not (Test-Path -LiteralPath $SetVars -PathType Leaf)) { throw "Intel oneAPI setvars not found: $SetVars; pass -SetVars or use -SkipBuild." }
    $script:lqmdExe = Join-Path $script:testRoot 'mock_lqmd.exe'
    $script:exfragExe = Join-Path $script:testRoot 'mock_exfrag.exe'
    Invoke-Ifx $Compiler $SetVars $script:lqmdExe $fixture
    Invoke-Ifx $Compiler $SetVars $script:exfragExe $fixture
} else {
    $script:lqmdExe = $LqmdMockExecutable
    $script:exfragExe = $ExfragMockExecutable
    foreach ($exe in @($script:lqmdExe,$script:exfragExe)) {
        if (-not (Test-Path -LiteralPath $exe -PathType Leaf)) { throw "With -SkipBuild, pass both compiled mock executables; missing: $exe" }
    }
}

$source = Join-Path $script:testRoot 'source'
$input = Join-Path $source 'Input file'
New-Item -ItemType Directory -Path $input, (Join-Path $source 'gdat'), (Join-Path $source 'tbl') -Force | Out-Null
Write-Utf8 (Join-Path $input 'LQMD.in') "1, 7, 0, 0 ! event range`r`n0, 0, 0 ! SURF_COAL`r`n8.25, 8.25, 0.5, 8.25, 1 ! x00 controls`r`n"
Write-Utf8 (Join-Path $input 'm_select.in') "0 ! m_select`r`n1 ! mass_framework`r`n314159 ! comparison_seed`r`n"
Write-Utf8 (Join-Path $input 'anahyp.in') @'
%%CONT_ANA
0 ! NTH_RUN
0 ! NTH_NOH
1, 7 ! NRMIN, NRMAX
%%Simulation
7 !events
8.25 !X00CUT
1, 1, 1 ! kktot, kkmin, kkmax
%%Models
1, 1 ! COLM, GEMINI
1, 10, 0 ! ANGM_DEC, ANGM_LIGHT, TESTANG
0 ! SURF_COAL
%%Yields
1, 1, 0, 75 ! CNPCS, IPCS, LCPs, TEELECS, TECOLD
'@
Write-Utf8 (Join-Path $input 'CNPCS_in.in') @'
name Z N ! Software smoke fixture species
--------------
225Ac 89 136
225Ra 88 137
227Ac 89 138
223Ra 88 135
227Th 90 137
--------------
'@
Write-Utf8 (Join-Path $source 'pini.dat') "MOCK INPUT`r`n"
Write-Utf8 (Join-Path $source 'tini.dat') "MOCK INPUT`r`n"
$sourceFilesBefore = @(Get-ChildItem -LiteralPath $source -File -Recurse | Sort-Object FullName | ForEach-Object { [pscustomobject]@{ Path=$_.FullName.Substring($source.Length); Hash=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash } })

$configPath = Join-Path $script:testRoot 'controller.json'
$mainRun = Join-Path $script:testRoot 'run main spaced'
$config = New-RunConfig $configPath $mainRun
$configArgs = @('-Stage','L','-Config',$configPath)
$controllerArgs = @('-File',$backend) + $configArgs
$prepareArgs = @('-File',$backend) + $configArgs + @('-PrepareOnly')
$prepareFirst = Invoke-Controller $prepareArgs
Assert-True ($prepareFirst.ExitCode -eq 0) "first PrepareOnly call succeeds: $($prepareFirst.Text)"
$prepareSecond = Invoke-Controller $prepareArgs
Assert-True ($prepareSecond.ExitCode -eq 0) "repeated PrepareOnly call succeeds: $($prepareSecond.Text)"
for ($i=1; $i -le 3; $i++) {
    $label = 'shard-{0:d4}' -f $i
    $attempts = @(Get-ChildItem -LiteralPath (Join-Path $mainRun "kinematics/$label") -Directory -Filter 'attempt-*')
    Assert-True ($attempts.Count -eq 1) "repeated PrepareOnly reuses LQMD shard $i attempt"
    $prepared = Get-Content -LiteralPath (Join-Path $attempts[0].FullName 'state.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    Assert-True ($prepared.status -eq 'prepared') "PrepareOnly leaves LQMD shard $i prepared"
}
$lResult = Invoke-Controller $controllerArgs
Assert-True ($lResult.ExitCode -eq 0) "LQMD mock stage failed: $($lResult.Text)"
$successfulPrepare = Invoke-Controller $prepareArgs
Assert-True ($successfulPrepare.ExitCode -eq 0) "PrepareOnly after success succeeds: $($successfulPrepare.Text)"
for ($i=1; $i -le 3; $i++) {
    $label = 'shard-{0:d4}' -f $i
    $attempts = @(Get-ChildItem -LiteralPath (Join-Path $mainRun "kinematics/$label") -Directory -Filter 'attempt-*')
    Assert-True ($attempts.Count -eq 1) "PrepareOnly after success does not create another LQMD shard $i attempt"
    $completed = Get-Content -LiteralPath (Join-Path $attempts[0].FullName 'state.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    Assert-True ($completed.status -eq 'completed') "PrepareOnly after success reuses completed LQMD shard $i"
}

$manifestPath = Join-Path $mainRun 'manifest.json'
Assert-True (Test-Path -LiteralPath $manifestPath) 'manifest must exist after LQMD stage'
$manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
$planned = @($manifest.shard_plan | Sort-Object index)
Assert-True ($planned.Count -eq 3) 'seven events split into exactly three shards'
$expectedRanges = @(@(1,3),@(4,5),@(6,7))
$expectedEvents = @(3,2,2)
$expectedSeeds = @()
$globalEvents = New-Object 'System.Collections.Generic.List[int]'
for ($i=0; $i -lt $planned.Count; $i++) {
    $shard = $planned[$i]
    Assert-True ([int]$shard.first -eq $expectedRanges[$i][0] -and [int]$shard.last -eq $expectedRanges[$i][1]) "shard $($i+1) has a contiguous non-overlapping global interval"
    Assert-True ([int]$shard.events -eq $expectedEvents[$i]) "shard $($i+1) has the exact tail event count"
    Assert-True ($expectedSeeds -notcontains [int]$shard.seed) "shard $($i+1) uses a distinct seed"
    $expectedSeeds += [int]$shard.seed
    $label = 'shard-{0:d4}' -f ($i + 1)
    $rawColl = Join-Path $mainRun ("kinematics/$label/attempt-0001/coll.dat")
    $stateFile = Join-Path $mainRun ("kinematics/$label/attempt-0001/state.json")
    $state = Get-Content -LiteralPath $stateFile -Raw -Encoding UTF8 | ConvertFrom-Json
    Assert-True ($state.status -eq 'completed') "LQMD shard $($i+1) completed"
    foreach ($line in [IO.File]::ReadAllLines($rawColl)) {
        if ($line -match '^mock-event\s+(\d+)$') { $globalEvents.Add([int]$Matches[1]) }
    }
}
Assert-True (($globalEvents | Sort-Object -Unique).Count -eq 7) 'all seven global event IDs occur once across shards'
Assert-True (($globalEvents | Sort-Object) -join ',' -eq '1,2,3,4,5,6,7') 'all global event IDs are present'

$aResult = Invoke-Controller (@('-File',$backend,'-Stage','A','-Config',$configPath))
Assert-True ($aResult.ExitCode -eq 0) "analysis mock stage failed: $($aResult.Text)"
foreach ($model in 1..5) {
    $merged = Join-Path $mainRun ("merged/model-$model/CNPCS_COLD.dat")
    Assert-True (Test-Path -LiteralPath $merged) "model $model merged CNPCS_COLD exists"
    $expected = ((3.0 * ($model * 10 + 1 * 5)) + (2.0 * ($model * 10 + 2 * 5)) + (2.0 * ($model * 10 + 3 * 5))) / 7.0
    $actual = Read-CnpcsValue $merged 'Z' 1
    Assert-Near $actual $expected 0.000001 "model $model merged yield is event-weighted across a short tail shard"
    Assert-True ([math]::Abs($actual - (($model * 10 + 5) + ($model * 10 + 10) + ($model * 10 + 15)) / 3.0) -gt 0.1) "model $model merge is not an unweighted shard mean"
    $lateOnlyExpected = 2.0 * ($model * 10 + 2 * 5) * 3.0 / 7.0
    $lateOnlyActual = Read-CnpcsValue $merged 'Z' 2
    Assert-Near $lateOnlyActual $lateOnlyExpected 0.000001 "model $model preserves a bin produced only by shard 2"
    foreach ($shardIndex in 1..3) {
        $label = 'shard-{0:d4}' -f $shardIndex
        $stateFile = Join-Path $mainRun ("analysis/model-$model/$label/attempt-0001/state.json")
        $state = Get-Content -LiteralPath $stateFile -Raw -Encoding UTF8 | ConvertFrom-Json
        Assert-True ($state.status -eq 'completed') "model $model analysis shard $shardIndex completed"
        Assert-True ([int]$state.seed -eq [int]$planned[$shardIndex - 1].seed) "model $model uses the same seed as paired shard $shardIndex"
    }
}
$targetRawOutput = Join-Path $mainRun 'analysis/model-1/shard-0001/attempt-0001/CNPCS_COLD.dat'
$originalRawOutput = [IO.File]::ReadAllBytes($targetRawOutput)
Clear-MergedOutputs $mainRun
Write-Utf8 $targetRawOutput "CORRUPTED AFTER WORKER COMPLETION`r`n"
$corruptText = Assert-ControllerRejected (@('-File',$backend,'-Stage','A','-Config',$configPath,'-MergeOnly')) 'a modified completed CNPCS output'
Assert-True ($corruptText -match 'Completed output was modified') 'the merge rejects a completed output whose digest changed'
Assert-True (-not (Test-Path -LiteralPath (Join-Path $mainRun 'merged/model-1/CNPCS_COLD.dat'))) 'the failed merge does not publish a corrupted CNPCS file'
[IO.File]::WriteAllBytes($targetRawOutput, $originalRawOutput)
Clear-MergedOutputs $mainRun
$null = Invoke-Controller (@('-File',$backend,'-Stage','A','-Config',$configPath,'-MergeOnly'))
Assert-True (Test-Path -LiteralPath (Join-Path $mainRun 'merged/model-5/CNPCS_HOT.dat')) 'MergeOnly can publish all previously completed analysis results'
$missingRawOutput = Join-Path $mainRun 'analysis/model-2/shard-0002/attempt-0001/CNPCS_HOT.dat'
[IO.File]::Delete($missingRawOutput)
Clear-MergedOutputs $mainRun
$missingText = Assert-ControllerRejected (@('-File',$backend,'-Stage','A','-Config',$configPath,'-MergeOnly')) 'a missing completed CNPCS output'
Assert-True (($missingText -match 'Completed output is missing') -or ($missingText -match 'missing')) 'the merge rejects a missing completed output'
Assert-True (-not (Test-Path -LiteralPath (Join-Path $mainRun 'merged/model-1/CNPCS_COLD.dat'))) 'the failed merge does not publish partial outputs'
$sourceFilesAfter = @(Get-ChildItem -LiteralPath $source -File -Recurse | Sort-Object FullName | ForEach-Object { [pscustomobject]@{ Path=$_.FullName.Substring($source.Length); Hash=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash } })
Assert-True (($sourceFilesBefore | ConvertTo-Json -Compress) -ceq ($sourceFilesAfter | ConvertTo-Json -Compress)) 'source input files are unchanged by all worker stages'

$badWorkerConfigPath = Join-Path $script:testRoot 'workers-nine.json'
$badWorkerConfig = New-RunConfig $badWorkerConfigPath (Join-Path $script:testRoot 'run-workers-nine') 7 3 9
$limitText = Assert-ControllerRejected (@('-File',$backend,'-Stage','L','-Config',$badWorkerConfigPath)) 'MaxWorkers=9'
Assert-True ($limitText -match 'safe range 1\.\.8') 'MaxWorkers=9 rejection states the 8-worker safety limit'

$changedConfig = New-RunConfig $configPath $mainRun 8 3 2
$changedText = Assert-ControllerRejected (@('-File',$backend,'-Stage','L','-Config',$configPath)) 'a changed event count for an existing run directory'
Assert-True ($changedText -match 'configuration or source hashes changed') 'changed configuration is rejected for an existing run directory'

$incompleteConfigPath = Join-Path $script:testRoot 'incomplete.json'
$incompleteRun = Join-Path $script:testRoot 'run-incomplete'
$null = New-RunConfig $incompleteConfigPath $incompleteRun
$null = Assert-ControllerRejected (@('-File',$backend,'-Stage','A','-Config',$incompleteConfigPath,'-MergeOnly')) 'merge before LQMD and analysis shards complete'

if (-not $SkipConcurrency) {
    $script:envSleepBefore = $env:MOCK_SLEEP_MS
    $env:MOCK_SLEEP_MS = '2500'
    try {
        $lockConfigPath = Join-Path $script:testRoot 'lock-first.json'
        $lockRun = Join-Path $script:testRoot 'run-lock-first'
        $null = New-RunConfig $lockConfigPath $lockRun 2 1 1
        $powershell = (Get-Command powershell.exe -ErrorAction Stop).Source
        $lockOut = Join-Path $script:testRoot 'lock-first.stdout.log'
        $lockErr = Join-Path $script:testRoot 'lock-first.stderr.log'
        $lockArgs = @('-NoProfile','-ExecutionPolicy','Bypass','-File',(Quote-Argument $backend),'-Stage','L','-Config',(Quote-Argument $lockConfigPath)) -join ' '
        $firstController = Start-Process -FilePath $powershell -ArgumentList $lockArgs -PassThru -WindowStyle Hidden -RedirectStandardOutput $lockOut -RedirectStandardError $lockErr
        $null = $firstController.Handle
        $startedFile = Join-Path $lockRun 'kinematics/shard-0001/attempt-0001/mock-started.txt'
        $deadline = [DateTime]::UtcNow.AddSeconds(15)
        while (-not (Test-Path -LiteralPath $startedFile)) {
            $firstController.Refresh()
            if ($firstController.HasExited -or [DateTime]::UtcNow -gt $deadline) { throw 'The first controller did not enter its mock worker.' }
            Start-Sleep -Milliseconds 50
        }
        $lockConfigSecond = Join-Path $script:testRoot 'lock-second.json'
        $lockRunSecond = Join-Path $script:testRoot 'run-lock-second'
        $null = New-RunConfig $lockConfigSecond $lockRunSecond 2 1 1
        $lockText = Assert-ControllerRejected (@('-File',$backend,'-Stage','L','-Config',$lockConfigSecond)) 'a concurrent independent controller'
        Assert-True ($lockText -match 'Another Windows batch controller is active') 'a concurrent controller is rejected by the global mutex'
        $firstController.WaitForExit()
        $firstController.Refresh()
        if ($firstController.ExitCode -ne 0) { throw "First controller failed during mutex test: $([IO.File]::ReadAllText($lockErr))" }
    } finally {
        if ($null -eq $script:envSleepBefore) { Remove-Item Env:MOCK_SLEEP_MS -ErrorAction SilentlyContinue }
        else { $env:MOCK_SLEEP_MS = $script:envSleepBefore }
    }
}

Write-Output "PASS: Windows controller software-only smoke tests completed under $script:testRoot"
Write-Output 'PASS: tail partition, unique event IDs/seeds, five-model paired seeds, weighted merge, source preservation and invalid-config rejection.'
Write-Output 'SCOPE: mock fixture outputs are software test data only; no physics or production 100000-event run was performed.'
