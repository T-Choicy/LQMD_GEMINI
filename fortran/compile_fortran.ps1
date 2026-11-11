param(
    [Parameter(Mandatory=$true)]
    [string]$SourceFile,
    [string]$OutputFile
)

# 1. 初始化 Intel oneAPI + Visual Studio 环境
Write-Host "Initializing Intel oneAPI environment..."
$setvarsLog = cmd /c "call ""`"C:\Program Files (x86)\Intel\oneAPI\setvars.bat`"" 2>&1 && echo SETVARS_OK"
if ($LASTEXITCODE -ne 0 -or $setvarsLog -notmatch "SETVARS_OK") {
    Write-Warning "setvars.bat had issues, trying manual PATH setup..."
    $ifxBin = "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin"
    if (Test-Path "$ifxBin\ifx.exe") {
        $env:PATH = "$ifxBin;$env:PATH"
        Write-Host "ifx PATH added manually."
    }
}

# 2. 找 ifx
$ifx = Get-Command ifx -ErrorAction SilentlyContinue
if (-not $ifx) {
    Write-Error "ifx not found! OneAPI environment setup failed."
    exit 1
}
Write-Host "Using: $($ifx.Source)"

# 3. 编译
if (-not $OutputFile) {
    $OutputFile = [System.IO.Path]::ChangeExtension($SourceFile, ".exe")
}
Write-Host "Compiling: $SourceFile -> $OutputFile"
$compileResult = & ifx $SourceFile -o $OutputFile 2>&1
if ($LASTEXITCODE -eq 0) {
    Write-Host "Compile SUCCESS: $OutputFile"
} else {
    Write-Error "Compile FAILED:"
    Write-Error $compileResult
    exit 1
}
