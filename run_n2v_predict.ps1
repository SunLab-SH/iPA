
conda activate ipa_env
cd "d:\Gitspace\ipa_full\iPA"
$env:PYTHONPATH = "d:\Gitspace\ipa_full\iPA"
python d:\Gitspace\ipa_full\iPA\examples\examples_sxt\demo_n2v_predict.py --model_path "d:\Gitspace\ipa_full\iPA\ipa\processing\denoising\models\n2v_3D\sxt_repro_run_final.pth" --n_filters 32
