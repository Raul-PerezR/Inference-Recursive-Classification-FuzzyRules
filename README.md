# Inference-Recursive-Classification-FuzzyRules
Software used at work "Use of search techniques to define efficient inference models for classification problems with a high number of fuzzy rules"

## How to compile?

cd src
make

## Sintax program?

Sintax:

 _InferStudy -e <path/seed problem> -model <num> [-nlabel <num>] [-sd <num>] [-d <num>] [-maxrules <num>] [-PerCentOnTest <real_num>]_
  
Parameters: 
	 -e  <path/seed problem> directory path of the problem and seed of the files 
	 -model <num> to select the inference model:
		1	Standard Inference
		2	Standard Inference Prunned
		3	Neighborhood Inference
		4	Heuristic Neighborhood Inference
		5	Heuristic Nearby Neighborhood Inference
		6	Hybrid Inference
	 -nlabel <num> number of labels used by discretize continuous variable. By default nlabel = 2  
	 -sd <num> seed for the random number generator. By default sd = 0 
	 -d <num> when model 5 or 6 is selected, this parameter establishes the maximum distance with the center rule. By default d = 0 
	 -maxrules <num> when model 3, 4, 5 or 6 is selected, this parameter fixes the limit in the number of rule for explored. By default maxrules = 1024
	 -PerCentOnTest <real_num> establishes de percentage of the examples from the test set on which the inference is applied. By default PerCentOnTest = 1.0
