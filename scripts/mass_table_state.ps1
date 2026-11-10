# Shared JSON state must permit Windows rename/delete while a reader is open.
function Read-MassRunJson {
    param([Parameter(Mandatory = $true)][string]$Path)
    for ($attempt = 0; $attempt -lt 40; $attempt++) {
        $stream = $null
        $reader = $null
        try {
            $share = [IO.FileShare]::ReadWrite -bor [IO.FileShare]::Delete
            $stream = [IO.File]::Open($Path, [IO.FileMode]::Open, [IO.FileAccess]::Read, $share)
            $reader = [IO.StreamReader]::new($stream, [Text.Encoding]::UTF8)
            return ($reader.ReadToEnd() | ConvertFrom-Json -ErrorAction Stop)
        } catch [IO.IOException] {
            if ($attempt -eq 39) { throw }
            Start-Sleep -Milliseconds 50
        } catch [UnauthorizedAccessException] {
            if ($attempt -eq 39) { throw }
            Start-Sleep -Milliseconds 50
        } finally {
            if ($null -ne $reader) { $reader.Dispose() }
            elseif ($null -ne $stream) { $stream.Dispose() }
        }
    }
}

function Write-MassRunJson {
    param([Parameter(Mandatory = $true)][string]$Path,
          [Parameter(Mandatory = $true)][object]$Value)
    $temporary = "$Path.$PID.tmp"
    [IO.File]::WriteAllText($temporary, ($Value | ConvertTo-Json -Depth 12), [Text.UTF8Encoding]::new($false))
    for ($attempt = 0; $attempt -lt 40; $attempt++) {
        try {
            [IO.File]::Move($temporary, $Path, $true)
            return
        } catch [IO.IOException] {
            if ($attempt -eq 39) { throw }
            Start-Sleep -Milliseconds 50
        } catch [UnauthorizedAccessException] {
            if ($attempt -eq 39) { throw }
            Start-Sleep -Milliseconds 50
        }
    }
}
