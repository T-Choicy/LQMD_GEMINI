# Result export for the NTH_RUN=0 Windows shard controller (PowerShell 5.1).
function Read-WindowsCnpcs {
    param([string]$Path)
    $rows = @{ Z=@{}; N=@{}; A=@{}; I=@{} }
    $section = ''
    foreach ($line in [IO.File]::ReadAllLines($Path)) {
        if ($line -match 'AS\s+A\s+FUNCTION\s+OF\s+([ZNA])') { $section=$Matches[1]; continue }
        if ($line -match 'CROSS\s+SECTION\s+OF\s+ISOTOPICS') { $section='I'; continue }
        if ($line -match '^\s*(\d+)\s+(\S+)\s*$') {
            if (-not $section) { throw "CNPCS data without section: $Path" }
            $id = [int]$Matches[1]
            $value = [double]::Parse($Matches[2].Replace('D','E').Replace('d','e'),[Globalization.CultureInfo]::InvariantCulture)
            if ([double]::IsNaN($value) -or [double]::IsInfinity($value) -or $value -lt 0) { throw "Invalid CNPCS value: $Path" }
            if ($rows[$section].ContainsKey($id)) { throw "Duplicate CNPCS bin: $Path $section $id" }
            $rows[$section][$id] = $value
        } elseif ($line.Trim() -and $line -notmatch 'NUCLIDE\s+CROSS|CROSS\s+SECTION|^-+$') {
            throw "Unrecognized CNPCS record: $Path : $line"
        }
    }
    foreach ($axis in @('Z','N','A')) { if (-not $rows[$axis].Count) { throw "Missing CNPCS $axis distribution: $Path" } }
    return $rows
}

function Get-WindowsSpecies {
    param([string]$Path)
    $rows = New-Object 'System.Collections.Generic.List[object]'
    $lines = [IO.File]::ReadAllLines($Path)
    for ($i=2; $i -lt $lines.Count; $i++) {
        $data = ($lines[$i] -split '!',2)[0].Trim()
        if (-not $data -or $data.StartsWith('---')) { break }
        $fields = $data -split '\s+'
        if ($fields.Count -ne 3) { throw "Invalid CNPCS_in.in record: $Path" }
        $rows.Add([pscustomobject]@{ Index=$rows.Count+1; Name=$fields[0]; Z=[int]$fields[1]; N=[int]$fields[2] })
    }
    if (-not $rows.Count) { throw 'CNPCS_in.in contains no species.' }
    return $rows.ToArray()
}

function Publish-WindowsCsv {
    param([string]$Path, [object[]]$Rows)
    $temp = $Path + '.tmp-' + [guid]::NewGuid().ToString('N')
    [IO.File]::WriteAllLines($temp,[string[]]@($Rows | ConvertTo-Csv -NoTypeInformation),$script:Utf8NoBom)
    if (Test-Path -LiteralPath $Path) {
        if ((Get-Sha256 $temp) -cne (Get-Sha256 $Path)) { throw "Refusing to overwrite modified result CSV: $Path" }
        [IO.File]::Delete($temp)
    } else { [IO.File]::Move($temp,$Path) }
}

function Export-WindowsControllerResults {
    param($Cfg)
    $distributions = New-Object 'System.Collections.Generic.List[object]'
    $targets = New-Object 'System.Collections.Generic.List[object]'
    foreach ($model in $Cfg.Models) {
        $base = Join-Path $Cfg.RunDirectory ('merged/model-{0}' -f $model)
        $species = @(Get-WindowsSpecies (Join-Path $base 'Input file/CNPCS_in.in'))
        foreach ($phase in @('HOT','COLD')) {
            $data = Read-WindowsCnpcs (Join-Path $base ('CNPCS_{0}.dat' -f $phase))
            foreach ($axis in @('A','Z','N')) {
                $maximum = if ($axis -eq 'A') { 300 } else { 150 }
                for ($bin=1; $bin -le $maximum; $bin++) {
                    $value = 0.0
                    if ($data[$axis].ContainsKey($bin)) { $value=$data[$axis][$bin] }
                    $distributions.Add([pscustomobject][ordered]@{
                        model=$model; model_name=$script:Models[[string]$model]; stage=$phase; variable=$axis; bin=$bin
                        cross_section_mb=$value.ToString('G17',[Globalization.CultureInfo]::InvariantCulture); events=$Cfg.Events
                    })
                }
            }
            if ($phase -eq 'COLD') {
                foreach ($target in $script:targets) {
                    $match = @($species | Where-Object { $_.Z -eq $target.Z -and $_.N -eq $target.N })
                    if ($match.Count -ne 1) { throw "Requested isotope $($target.Name) must occur exactly once in CNPCS_in.in." }
                    $index = [int]$match[0].Index
                    if (-not $data.I.ContainsKey($index)) { throw "Missing isotope index $index in model $model COLD output." }
                    $targets.Add([pscustomobject][ordered]@{
                        model=$model; model_name=$script:Models[[string]$model]; isotope=$target.Name; Z=$target.Z; N=$target.N; A=$target.A
                        cross_section_mb=$data.I[$index].ToString('G17',[Globalization.CultureInfo]::InvariantCulture); events=$Cfg.Events
                    })
                }
            }
        }
    }
    Publish-WindowsCsv (Join-Path $Cfg.RunDirectory 'merged/distributions.csv') $distributions.ToArray()
    Publish-WindowsCsv (Join-Path $Cfg.RunDirectory 'merged/target_isotopes.csv') $targets.ToArray()
}
