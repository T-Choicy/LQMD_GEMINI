[CmdletBinding()]
param(
    [string]$Compiler,
    [string]$SetVars,
    [switch]$Rebuild
)

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$fortranRoot = Join-Path $repoRoot 'fortran'
$outputRoot = Join-Path $fortranRoot 'x64/Release/win_batch'
$taskTemp = Join-Path $repoRoot 'tmp/build/windows-controllers'

function Find-OneApiCompiler {
    $installRoots = @()
    foreach ($base in @(${env:ProgramFiles(x86)}, $env:ProgramFiles)) {
        if ($base) {
            $installRoots += Join-Path $base 'Intel/oneAPI'
        }
    }
    foreach ($installRoot in $installRoots) {
        $compilerRoot = Join-Path $installRoot 'compiler'
        if (-not (Test-Path -LiteralPath $compilerRoot -PathType Container)) {
            continue
        }
        foreach ($versionDir in (Get-ChildItem -LiteralPath $compilerRoot -Directory |
                    Sort-Object Name -Descending)) {
            $candidate = Join-Path $versionDir.FullName 'bin/ifx.exe'
            if (Test-Path -LiteralPath $candidate -PathType Leaf) {
                return [pscustomobject]@{
                    Compiler = $candidate
                    SetVars = Join-Path $installRoot 'setvars.bat'
                }
            }
        }
    }
    throw 'Intel oneAPI ifx.exe was not found. Supply -Compiler and optionally -SetVars.'
}

if (-not $Compiler) {
    $installation = Find-OneApiCompiler
    $Compiler = $installation.Compiler
    if (-not $SetVars) { $SetVars = $installation.SetVars }
}
$Compiler = [IO.Path]::GetFullPath($Compiler)
if (-not (Test-Path -LiteralPath $Compiler -PathType Leaf)) {
    throw "Fortran compiler not found: $Compiler"
}
if ($SetVars) {
    $SetVars = [IO.Path]::GetFullPath($SetVars)
    if (-not (Test-Path -LiteralPath $SetVars -PathType Leaf)) {
        throw "oneAPI environment script not found: $SetVars"
    }
}

New-Item -ItemType Directory -Path $taskTemp -Force | Out-Null
$targets = @(
    @{ Name = 'ctrlL_win'; Source = 'ctrlL_win.f' },
    @{ Name = 'ctrlA_win'; Source = 'ctrlA_win.f' },
    @{ Name = 'LQMD'; Source = 'LQMD.f' },
    @{ Name = 'exfrag'; Source = 'exfrag.f' }
)
$targetFiles = foreach ($target in $targets) {
    Join-Path $outputRoot ($target.Name + '.exe')
}
$outputRoot = [IO.Path]::GetFullPath($outputRoot)
$approvedOutput = [IO.Path]::GetFullPath((Join-Path $fortranRoot 'x64/Release/win_batch'))
if (-not [string]::Equals($outputRoot,$approvedOutput,[StringComparison]::OrdinalIgnoreCase)) {
    throw "Resolved build output is outside the dedicated win_batch directory: $outputRoot"
}
if (Test-Path -LiteralPath $outputRoot -PathType Container) {
    $existing = @($targetFiles | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf })
    if ($existing.Count -eq 4 -and -not $Rebuild) {
        Write-Host "All four win_batch executables already exist in $outputRoot."
        Write-Host 'Use -Rebuild to archive them and compile new binaries.'
        exit 0
    }
    if (-not $Rebuild) {
        throw 'Partial win_batch output exists. Use -Rebuild to archive it before compiling.'
    }
    $active = @(Get-Process -Name ctrlL_win,ctrlA_win,LQMD,exfrag,Console1 -ErrorAction SilentlyContinue)
    if ($active.Count -gt 0) {
        throw 'Cannot rebuild while a controller or physics executable is running.'
    }
    $archiveRoot = Join-Path $taskTemp ('win_batch-archive-' +
        [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffZ') + '-' +
        [guid]::NewGuid().ToString('N'))
    if (Test-Path -LiteralPath $archiveRoot) { throw "Archive already exists: $archiveRoot" }
    Move-Item -LiteralPath $outputRoot -Destination $archiveRoot
    Write-Host "Existing build outputs preserved at $archiveRoot"
}
New-Item -ItemType Directory -Path $outputRoot -Force | Out-Null

function Quote-BatchArgument([string]$Value) {
    if ($Value.Contains('"') -or $Value.Contains('%')) {
        throw "Build path contains an unsupported quote or percent sign: $Value"
    }
    return '"' + $Value + '"'
}

$batchPath = Join-Path $taskTemp 'compile_windows_controllers.cmd'
$batchLines = [Collections.Generic.List[string]]::new()
$batchLines.Add('@echo off')
if ($SetVars) {
    $batchLines.Add('call ' + (Quote-BatchArgument $SetVars) + ' >nul')
    $batchLines.Add('if errorlevel 1 exit /b %errorlevel%')
}

foreach ($target in $targets) {
    $sourcePath = Join-Path $fortranRoot $target.Source
    $targetObj = Join-Path $outputRoot ('obj/' + $target.Name)
    $modulePath = Join-Path $targetObj 'mod'
    New-Item -ItemType Directory -Path $modulePath -Force | Out-Null
    $objectPath = Join-Path $targetObj ($target.Name + '.obj')
    $executablePath = Join-Path $outputRoot ($target.Name + '.exe')
    $args = @(
        '/nologo', '/O2', '/fixed', '/extend-source:132',
        '/check:bounds', '/traceback', '/libs:static',
        ('/I' + (Quote-BatchArgument $fortranRoot)),
        ('/module:' + (Quote-BatchArgument $modulePath)),
        ('/object:' + (Quote-BatchArgument $objectPath)),
        ('/exe:' + (Quote-BatchArgument $executablePath)),
        (Quote-BatchArgument $sourcePath)
    )
    $batchLines.Add('call ' + (Quote-BatchArgument $Compiler) + ' ' +
        ($args -join ' '))
    $batchLines.Add('if errorlevel 1 exit /b %errorlevel%')
}
$batchLines.Add('exit /b 0')
[IO.File]::WriteAllLines($batchPath, $batchLines, [Text.Encoding]::ASCII)

Write-Host "Compiler: $Compiler"
if ($SetVars) { Write-Host "Environment: $SetVars" }
Write-Host "Building four x64 programs in $outputRoot"
$argumentText = '/d /s /c "' + $batchPath + '"'
$process = Start-Process -FilePath $env:ComSpec -ArgumentList $argumentText `
    -WorkingDirectory $fortranRoot -PassThru -WindowStyle Hidden
$null = $process.Handle
# Wait for the compiler batch itself, without waiting on its descendant job.
$process.WaitForExit()
$process.Refresh()
if ($process.ExitCode -ne 0) {
    throw "ifx build failed with exit code $($process.ExitCode). See $batchPath."
}

foreach ($file in $targetFiles) {
    if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
        throw "Compiler returned success but did not create $file"
    }
}
Write-Host 'Build completed: ctrlL_win.exe, ctrlA_win.exe, LQMD.exe, exfrag.exe.'
