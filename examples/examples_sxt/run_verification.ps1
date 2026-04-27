
Write-Host "Running N2N Prediction using ipa_env..."
cmd /c "conda run -n ipa_env python demo_SXT_n2n_denoising.py > n2n_log.txt 2>&1"
if ($LASTEXITCODE -ne 0) { Write-Error "N2N Prediction Failed" }

Write-Host "Running N2V Prediction using ipa_env..."
cmd /c "conda run -n ipa_env python demo_SXT_denoiseing.py > n2v_log.txt 2>&1"
if ($LASTEXITCODE -ne 0) { Write-Error "N2V Prediction Failed" }

Write-Host "Verification Complete."
