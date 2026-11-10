[CmdletBinding()]
param([switch]$Probe,[switch]$Abandon)
$ErrorActionPreference='Stop'
$script:WindowsHost=[Environment]::OSVersion.Platform -eq [PlatformID]::Win32NT
if($script:WindowsHost){throw 'This test exercises the Linux file-lock branch; use the Windows concurrency suite on Windows.'}
$script:HeldSlots=@{}
$backend=Join-Path (Split-Path -Parent $PSScriptRoot) 'scripts/windows_event_controller.ps1'
$ast=[Management.Automation.Language.Parser]::ParseFile($backend,[ref]$null,[ref]$null)
foreach($name in @('Get-ComputeSlot','Close-ComputeSlot')){
 $function=$ast.FindAll({param($node) $node -is [Management.Automation.Language.FunctionDefinitionAst]},$false)|Where-Object Name -eq $name
 if(@($function).Count -ne 1){throw "Missing production function: $name"}
 Invoke-Expression $function.Extent.Text
}
if($Probe){
 $slot=Get-ComputeSlot -NoWait
 if($slot){Write-Output 'acquired';if(-not $Abandon){Close-ComputeSlot $slot}}else{Write-Output 'blocked'}
 exit 0
}
$slots=[Collections.Generic.List[object]]::new()
try{
 foreach($i in 1..25){$slot=Get-ComputeSlot -NoWait;if(-not $slot){throw "Missing free slot $i"};$slots.Add($slot)}
 if(Get-ComputeSlot -NoWait){throw 'The same process acquired more than 25 slots.'}
 if((& (Join-Path $PSHOME 'pwsh') -NoProfile -File $PSCommandPath -Probe) -ne 'blocked'){throw 'Another process bypassed the 25-slot cap.'}
 Close-ComputeSlot $slots[0];$slots.RemoveAt(0)
 if((& (Join-Path $PSHOME 'pwsh') -NoProfile -File $PSCommandPath -Probe -Abandon) -ne 'acquired'){throw 'Released slot unavailable to another process.'}
 $slot=Get-ComputeSlot -NoWait
 if(-not $slot){throw 'Process exit failed to release its file lock.'}
 $slots.Add($slot)
 Write-Output 'PASS: 25 distinct local slots, cross-process exclusion, explicit release and process-exit release.'
}finally{foreach($slot in $slots){Close-ComputeSlot $slot}}
