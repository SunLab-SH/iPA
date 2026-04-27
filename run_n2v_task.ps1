
conda activate ipa_env
cd "d:\Gitspace\ipa_full\iPA"
$env:PYTHONPATH = "d:\Gitspace\ipa_full\iPA"
python -m ipa.processing.denoising.noise2void.run_train --input "d:\Gitspace\ipa_full\iPA\data\sxt_images\Stevens_pancreatic_INS_1E_784_5_pre_rec.mrc" --name sxt_repro_run --epochs 5 --n_filters 32 --batch_size 2
