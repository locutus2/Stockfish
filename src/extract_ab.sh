#!/bin/bash
sed 's/balanced accuracy/balanced_accuracy/; s/false positive rate/false_negative_rate/; s/false negative rate/false_positive_rate/' \
	| gawk ' \
	BEGIN { count = 0; } \
	$0 ~ /Iteration/ { count++; } \
	$2 ~ /support:/ { sup[count] = $3; } \
	$2 ~ /^accuracy/ { acc[count] = $3; } \
	$2 ~ /balanced/ { balacc[count] = $3; } \
	$2 ~ /true_positive/ { pos[count] = $3; } \
	$2 ~ /true_negative/ { neg[count] = $3; } \
	$2 ~ /false_positive/ { pos[count] = (100-($3 + 0)) "%"; } \
	$2 ~ /false_negative/ { neg[count] = (100-($3 + 0)) "%"; } \
	$2 ~ /precision/ { prec[count] = $3; } \
	$2 ~ /F1/ { F1[count] = $3; } \
	$2 ~ /MCC/ { MCC[count] = $3; } \
	END { print "Iteration Support Accuracy BalancedAccuracy TruePositiveRate TrueNegativeRate Precision F1-Score MCC"; \
		for(c in acc) {\
		print c, sup[c], acc[c], balacc[c], pos[c], neg[c], prec[c], F1[c], MCC[c]; \
		} } \
     ' | tr '. ' ',;'
