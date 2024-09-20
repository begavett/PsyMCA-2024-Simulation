9/20/2024

Brandon’s interim work


Another discretization scheme: Doug Tommet and Rich Jones are, for this and other projects, converting a Stata function to R: Paul Crane’s function with Gibbons and Jones, which interfaced with parscale IRT 


Something runs and works!
 Simulate continuous variables, discretize according to schemes, and run with IRT to determine whether DIF is artificially induced

Instead of using real scales, use a synthetic scale and convert to more typical neurospych data. 
Only real data: HAALSI for mean and SD memory score factor groups to set population parameters

See Brandon’s code and annotations for specifications

Dialogue about tradeoffs between weak categorical scales from real data and a continuous synthetic scale

Dan: should be tweaked to reflect severe skew usually observed in real data, manipulating distributions to evaluate impact

Doug: instead of a continuous z-score parameter, could instead simulate categorical variables

Doug: how would you discretize if there’s little overlap between samples? Based on the higher or lower scoring sample? 

Need to replicate issues of overlap seen in real studies without introducing DIF

Universal lookup tables

Equal interval, equal frequency, kmeans, manual, and collapse


DIF detection using continuous variables as a benchmark as Jones has done: mimic modeling. Then can better quantify false positives

Brandon has more code to share

Actionable tasks to divvy up work: manual recoding down the road and data distributions, work on coding including mimic modeling or MI testing, continuous item operationalization, need to make data resemble actual data

Before next week: look at the code!
 https://github.com/begavett/PsyMCA-2024-Simulation/blob/main/simulate-dif-detection.R


9/27 
Naive Kendra question: could any parameter, data distribution, etc. come from HCAP studies to ground the implications of this project in real data? 
