# Inference-Recursive-Classification-FuzzyRules
Software used at work "Use of search techniques to define efficient inference models for classification problems with a high number of fuzzy rules"

## How to compile?

cd src

make

## How to run?

Sintax:

`InferStudy -e <path/seed problem> -model <num> [-nlabel <num>] [-sd <num>] [-d <num>] [-maxrules <num>] [-PerCentOnTest <real_num>]`
  
Parameters: 
* `-e  <path/seed problem>` directory path of the problem and seed of the files 
* `-model <num>` to select the inference model:
  * `1`	Standard Inference
  * `2`	Standard Inference Prunned
  * `3`	Neighborhood Inference
  * `4`	Heuristic Neighborhood Inference
  * `5`	Heuristic Nearby Neighborhood Inference
  * `6`	Hybrid Inference
* `-nlabel <num>` number of labels used by discretize continuous variable. By default nlabel = 2  
* `-sd <num>` seed for the random number generator. By default sd = 0 
* `-d <num>` when model 5 or 6 is selected, this parameter establishes the maximum distance with the center rule. By default d = 0 
* `-maxrules <num>` when model 3, 4, 5 or 6 is selected, this parameter fixes the limit in the number of rule for explored. By default maxrules = 1024
* `-PerCentOnTest <real_num>` establishes de percentage of the examples from the test set on which the inference is applied. By default PerCentOnTest = 1.0


## Some examples

##### ./InferStudy -e ../databases/texture/texture -model 1 -nlabel 5
Result of 10 cross-validation on the *texture* database using the *Standard Inference* and 5 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables.

##### ./InferStudy -e ../databases/census/census -model 4 -nlabel 5 -maxrules 1024
Result of 10 cross-validation on the *census* database using the *Neighborhood Inference Heuristic* and 5 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables and setting the parameter on the maximum number of rule to search in the example neighbor to 2^10 = 1024.


##### ./InferStudy -e ../databases/sonar/sonar -model 6 -nlabel 3 -maxrules 1024 -d 2
Result of 10 cross-validation on the *sonar* database using the *Hybrid Inference* and 3 uniformly distributed labels cutoff at 0.5 over the universe of discourse of continuous variables and setting the parameter on the maximum number of rule to search in the example neighbor to 2^10 = 1024 and setting the parameter on the distance on the center rules to 2.

## Using script *Launch_all.sh* for replicating the experimentation
The *kddcup*, *susy*, *higgs* and *hepmass* databases have not been included in this repository due to their excessive size. You can contact us to obtain these databases just as they were used for this experimentation.

### Obtaining Table 2: Results on training sets
The results shown for the training phase are common to all the models studied. For that reason, it is the same on which inference model is run and the process is faster if inference is not done on the test set. Therefore, we set the parameter `-PerCentOnTest 0 `. Thus, the results shown in Table 2 can be obtained by running the script setting as parameters after each database `-model 1 -nlabel 5 -PerCentOnTest 0`. 

By example, for *census* database it will be `./InferStudy -e ../databases/census/census -model 4 -nlabel 5 -PerCentOnTest 0`.


### Obtaining results of Algorithm 2 for Tables 3,4,5,6,7 and 8
`-model 2 -nlabel 5` are the parameters for all databases except by *higgs* and *hepmass*. In the last two cases we must to include the PerCentOnTest parameter with 0.001 value.

By example, for *census* database (and the rest of databases except *higgs* and *hepmass*), we will use `./InferStudy -e ../databases/census/census -model 2 -nlabel 5`. In the case of *higgs* and *hepmass*, we will use  `./InferStudy -e ../databases/higgs/higgs -model 2 -nlabel 5 -PerCentOnTest 0.001` and `./InferStudy -e ../databases/higgs/higgs -model 2 -nlabel 5 -PerCentOnTest 0.001`.


### Obtaining results of Algorithm 3 and 4 for Tables 3,4,5,6,7 and 8
In the case of algorithms 3 and 4, the way to call them and modify the script is the same as described for algorithm 2, including the same exceptions on the *higgs* and *hepmass* databases. The differences are that in algorithm 3 the value 3 is associated to the parameter *model* and in algorithm 4 the value 4. The other important difference is that is need to fix a value of a `maxrules` parameter. In the case of this experimentation, it was fixed to `1024`

An example, for *census* database (and the rest of databases except *higgs* and *hepmass*) and algorithm 3, we will use `./InferStudy -e ../databases/census/census -model 3 -nlabel 5 -maxrules 1024`. In the case of *higgs* and *hepmass* is appened `-PerCentOnTest 0.001`.

### Obtaining Tables 9 and 10: Results on Algorithm 5 with different values of distance
Table 9 shows the results of a modification of Algorithm 4 where the space where the best rule is searched in the neighborhood of the example is restricted. The values 0,1,2,3,5 have been considered for this distance parameter. In this case, the script must be invoked 5 times, changing in each case the value associated to the parameter `-d` from 0 to 5 without considering the 4. Thus, for all the databases, and distance 0, this will be the way to invoke `-model 5 -nlabel 5 -d 0`.

In the same way, just replace the 0 of the parameter `-d` by a value in {1,2,3,5} to obtain the result with the changed distance.  To these parameters, in the case of *higgs* and *hepmas* it is necessary to add the parameter `-PerCentOnTest 0.001`.

Only 16 of the 50 databases are used in that part of the paper. You could comments in the script non-used databases for accelerating the process.


### Obtaining Table 11: Results on Algorithm 6

Similarly as described in the previous section, only 16 of the 50 databases have been used in this part. You could put in comment in the script those unused databases.

The way to invoke this version is a combination of the invocation between algorithms 4 and 5. For this reason, both the `-maxrules` and `-d` parameters are used. The general modification in the script for each database (and after setting with the `-e` parameter the database) is `-model 6 -nlabel 5 -maxrules 1024 -d 2`.As in all previous cases, for the *higgs* and *hepmass* databases `-PerCentOnTest 0.001` must be added.
