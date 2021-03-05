# Inference-Recursive-Classification-FuzzyRules
Software used at work "Use of search techniques to define efficient inference models for classification problems with a high number of fuzzy rules"

## How to compile?

cd src
make

## How to run?

Sintax:

##### InferStudy -e <path/seed problem> -model <num> [-nlabel <num>] [-sd <num>] [-d <num>] [-maxrules <num>] [-PerCentOnTest <real_num>]
  
Parameters: 
* -e  <path/seed problem> directory path of the problem and seed of the files 
* -model <num> to select the inference model:
  * 1	Standard Inference
  * 2	Standard Inference Prunned
  * 3	Neighborhood Inference
  * 4	Heuristic Neighborhood Inference
  * 5	Heuristic Nearby Neighborhood Inference
  * 6	Hybrid Inference
* -nlabel <num> number of labels used by discretize continuous variable. By default nlabel = 2  
* -sd <num> seed for the random number generator. By default sd = 0 
* -d <num> when model 5 or 6 is selected, this parameter establishes the maximum distance with the center rule. By default d = 0 
* -maxrules <num> when model 3, 4, 5 or 6 is selected, this parameter fixes the limit in the number of rule for explored. By default maxrules = 1024
* -PerCentOnTest <real_num> establishes de percentage of the examples from the test set on which the inference is applied. By default PerCentOnTest = 1.0


## Some examples

##### ./InferStudy -e ../databases/texture/texture -model 1 -nlabel 5
Result of 10 cross-validation on the *texture* database using the *Standard Inference* and 5 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables.

##### ./InferStudy -e ../databases/census/census -model 4 -nlabel 5 -maxrules 1024
Result of 10 cross-validation on the *census* database using the *Neighborhood Inference Heuristic* and 5 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables and setting the parameter on the maximum number of rule to search in the example neighbor to 2^10 = 1024.


##### ./InferStudy -e ../databases/sonar/sonar -model 6 -nlabel 3 -maxrules 1024 -d 2
Result of 10 cross-validation on the *sonar* database using the *Hybrid Inference* and 3 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables and setting the parameter on the maximum number of rule to search in the example neighbor to 2^10 = 1024 and setting the parameter on the distance on the center rules to 2.

## Using script *Launch_all.sh* for replicating the experimentation
The *kddcup*, *susy*, *higgs* and *hepmass* databases have not been included in this repository due to their excessive size. You can contact us to obtain these databases just as they were used for this experimentation.

### Table 2: Result on training sets
* The results shown for the training phase are common to all the models studied. For that reason, it is the same on which inference model is run and the process is faster if inference is not done on the test set. Therefore, we set the parameter *PerCentOnTest* to 0. Thus, the results shown in Table 2 can be obtained by running the script setting as parameters after each database -model 1 -nlabel 5 -PerCentOnTest 0.
