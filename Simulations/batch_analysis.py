"""
read csv
for each row:
	calculate B/T, paramD/paramB mix and match, parami-paramf, chi, model-truth_residuals,
save csv 
decide conditions for a good fit, disputed-fit, bad-fit
	delta_params, std_errors, redchi
	good=matches truth with low error, disputed=(low error but mismatched,)OR(high error and matched) bad=high/no errors and mismatched

	low_error = error<1%, good_match= p_fit-error<p_truth<p_fit+error

	good if:
		error < 1%
		good_match

	disputed if:
		error > 1% OR bad_match

	else: bad

classify all fits
"""

