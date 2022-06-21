# segmentation_shared
Public repo for Optimization of Portfolio Partitioning for Chain-Ladder Reserve Estimates article.
Prototype of chain-ladder segmentation optimizer
The main function runner optimizes segmentation ussing:
a) Greedy algorithm
b) Exhaustive algorithm


## Running 
Unpack zip file and its subdirectories in to a working directory.
Change the path of the working directory <workdir> at the beginning of the segmentation_runner_vXX.R and run the script.
Run the segmentation_runner_vXX.R script.

## Input
- All names and algorithm parameters are stored in config_file in <workdir> directory
- The application will process a <input_data_file> stored in the `<workdir>/<wd_data>` directory specified in congig.
- The application uses functions from the `<workdir>/<wd_functions>` directory. 
All these directories are set in <config_file>

## Output
- Output is displayed in console and as a interactive plot (for greedy algorithm only)
Example of an output:
#For both exhaustive as well a a greedy algorithm, a table of the optimal result is printed: 
                          rule  	  	               err      loss likelihood   penalty leaf_cnt
1 LoBcnt_1==1|LoBcnt_4==1 , LoBcnt_2==1 , LoBcnt_3==1 1.910843e+12 0.9191702  0.7542407 0.1649294        3

For greedy algorithm, a scheme is printed 
 Total              
  ¦--LoBcnt_2==0    
  ¦   ¦--LoBcnt_3==0
  ¦   °--LoBcnt_3==1
  °--LoBcnt_2==1 
and a graphical tree plot containing also relative volume of the nodes are added. 

Optimized partitioning
rlue ... partition rule ("," is the separator between segments). Syntax reflects one hot encoding syntax of the <categ_var_exhaustive> variable, e.g. LoBcnt_1==1 means LoBcnt == 1.
err ... error of the total estimated reserve evalueated on the future data (if contained in the dataset)
loss ... Value of the loss function (loss =likelihood + penalty)
leaf_cnt ... Number of segments in this partitioning.


##<input_data_file>
Input data file is stored as a csv and corresponds with usual reserving data input in a long format containing obligatory:
'id'...Claim id 
'occ_yr'...Occurrence year
'dev_yr'...Development year (starting at year 1)
<size_var>...Variable containing the triangle variable for which the reserve is to be calculated (paid, incurred amount etc...). Its name is set in <config_file>. 

Furthermore for greedy algorithm, there are three types of segmentation variables with user defined names. Assignment of the variables to each of these groups is done in <config_file>:
1)<metric_var>... Metric (numerical) variables. No encoding will be performed.
2) <categ_var_enc>...Categorical variables to be encoded using product of all estimated chin-ladder development factors as numerical values for each category. (Recommended if number of categories is higher than approximatelly 15.)
3) <categ_var>...Categorical variables for which all possible partitioning are created and set unique elements of this partitionings are encoded using one-hot-encoding (dummy variables). (Recommended if number of categories does not exceed approximatelly 15.)

For exhaustive search algorithm, there is only one categorical segmentation variable <categ_var_exhaustive>. If combinations of values of several categorical variables are to be used, combinations must be prepared outside the application.



##Configuration file <config_file>
Config file contains all names and algorithm parameters. 

<wd_data> ... Relative path to data (string)
<wd_functions> ... Relative path to application functions (string)
<input_data_file> ... Name of the input file (string)
<set_algorithm> ... A string value to switch between the optimization algorithms. Possibilities:"exhaustive" or "greedy". These algorithms require some different inputs.
<size_var> ... Name of the triangle variable for which the reserve is to be calculated.(string)
<categ_var_exhaustive> ... Only for exhaustive algorithm. A categorical variable, for which all potential partitionings are to be created and evalueated.(string)
<metric_var>, <categ_var_enc>, <categ_var> ...  String vectors reopresenting categorization of segmentation variables described in section "input_data_file" above. 
<use_set_var> ... Vector of strings defining subset of variables selected to be used for optimization. All variables selected must be contained in exactly one of the three categories <metric_var>, <categ_var_enc>, <categ_var>.

<trunc_penalty_mean> ... A numeric value between [0-1]. Use truncation of a given percentage when calculating mean penalty term. 
<penalty_fun> ... A string value. Set to "median" if median is desired instead of (truncated) mean in penalty term.  Use if mean penalty is numerically unstable. If not defined, truncated mean is used.
<variance_power> ... An integer used as the power parameter for variance function (Var = sigma * mu^variance_power). Recommended value (consistent with Mack's model) is 1. 
<ignore_dev> ...An integer defining how many of the corner values should be omitted in triangles. Recommendend value is 2.

<prob_grid> ... Numerical vector defining a grid of percentiles for numerical variables over which optimal value is searched for. Relevant for greedy algorithm only.
<min_volume> ... A stopping criterion of the optimization. Algorithm does not split if node volume is smaller than this value.
<calc_penalty> ... A logical value. Set "T" to include also penalty term.
<calc_future> ... A logical value. Set "T" to calculate also error for future predictions. If the value is "T", future data (=lower right triangle data) must be included in the input data.
