<#
.SYNOPSIS
Parse timing information out of LAMMPS stdout/log output.

.DESCRIPTION
LAMMPS prints, for every `run` command, a block that looks like:

    Loop time of 0.0367205 on 1 procs for 250 steps with 4000 atoms

    Performance: 11762228.140 tau/day, 6807.532 timesteps/s, 27.230 Matom-step/s
    99.7% CPU use with 1 MPI tasks x no OpenMP threads

    MPI task timing breakdown:
    Section |  min time  |  avg time  |  max time  |%varavg| %total
    ---------------------------------------------------------------
    Pair    | 0.023019   | 0.023019   | 0.023019   |   0.0 | 62.69
    Neigh   | 0.0091632  | 0.0091632  | 0.0091632  |   0.0 | 24.96
    Comm    | 0.0017677  | 0.0017677  | 0.0017677  |   0.0 |  4.81
    Output  | 7.2e-05    | 7.2e-05    | 7.2e-05    |   0.0 |  0.20
    Modify  | 0.0017543  | 0.0017543  | 0.0017543  |   0.0 |  4.78
    Other   |            | 0.0009463  |            |       |  2.58

ConvertFrom-LammpsLog turns this into a PSCustomObject so it can be
aggregated across repeated runs and exported as JSON/CSV.
#>

$script:SectionNames = @('Pair', 'Bond', 'Angle', 'Dihedral', 'Improper', 'Kspace',
    'Neigh', 'Comm', 'Output', 'Modify', 'Other')

function ConvertFrom-LammpsLog {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory, ValueFromPipeline)]
        [AllowEmptyString()]
        [string[]]$InputObject
    )
    begin {
        $lines = [System.Collections.Generic.List[string]]::new()
    }
    process {
        foreach ($line in $InputObject) { $lines.Add($line) }
    }
    end {
        $text = ($lines -join "`n")

        $result = [ordered]@{
            Version = $null
            Runs    = @()
        }

        if ($text -match '(?m)^LAMMPS \(([^)]*)\)') {
            $result.Version = $Matches[1].Trim()
        }

        $startMatches = [regex]::Matches($text, 'Loop time of')
        $blocks = @()
        for ($i = 0; $i -lt $startMatches.Count; $i++) {
            $start = $startMatches[$i].Index
            if ($i + 1 -lt $startMatches.Count) {
                $end = $startMatches[$i + 1].Index
            }
            else {
                $end = $text.Length
            }
            $blocks += $text.Substring($start, $end - $start)
        }

        foreach ($block in $blocks) {
            $run = [ordered]@{}

            if ($block -match 'Loop time of ([\d.eE+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms') {
                $run.LoopTime = [double]$Matches[1]
                $run.NProcs = [int]$Matches[2]
                $run.NSteps = [int]$Matches[3]
                $run.NAtoms = [int]$Matches[4]
            }

            if ($block -match '(?m)^Performance:\s*(.+)$') {
                $perf = $Matches[1].Trim()
                $run.Performance = $perf
                if ($perf -match '([\d.]+)\s*timesteps/s') {
                    $run.TimestepsPerSec = [double]$Matches[1]
                }
            }

            if ($block -match '([\d.]+)% CPU use with (\d+) MPI tasks? x (.+?) threads?') {
                $run.CpuPct = [double]$Matches[1]
                $run.MpiTasks = [int]$Matches[2]
                $threads = $Matches[3].Trim()
                if ($threads -match '^no') {
                    $run.OmpThreads = 0
                }
                elseif ($threads -match '^(\d+)') {
                    $run.OmpThreads = [int]$Matches[1]
                }
                else {
                    $run.OmpThreads = $null
                }
            }

            $sectionLines = $block -split "`n"
            $sections = [ordered]@{}
            for ($i = 0; $i -lt $sectionLines.Count; $i++) {
                if ($sectionLines[$i].Trim().StartsWith('Section |')) {
                    for ($j = $i + 2; $j -lt $sectionLines.Count; $j++) {
                        $row = $sectionLines[$j]
                        if ($row -notmatch '\|') { break }
                        $cells = $row -split '\|' | ForEach-Object { $_.Trim() }
                        if ($cells.Count -lt 6) { break }
                        $name = $cells[0]
                        if ($script:SectionNames -notcontains $name) { break }
                        $avg = $cells[2]
                        $pct = $cells[5]
                        $sections[$name] = [ordered]@{
                            AvgTime  = if ($avg) { [double]$avg } else { $null }
                            PctTotal = if ($pct) { [double]$pct } else { $null }
                        }
                    }
                    break
                }
            }
            $run.Sections = $sections

            $result.Runs += [PSCustomObject]$run
        }

        if ($result.Runs.Count -gt 0) {
            $last = $result.Runs[-1]
            foreach ($prop in $last.PSObject.Properties) {
                $result[$prop.Name] = $prop.Value
            }
        }

        return [PSCustomObject]$result
    }
}

Export-ModuleMember -Function ConvertFrom-LammpsLog
