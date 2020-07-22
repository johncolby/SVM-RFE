# (multiple) Support Vector Machine Recursive Feature Elimination (mSVM-RFE)

This package contains an [R](http://www.r-project.org) implementation of the mSVM-RFE algorithm ([Duan et al., 2005](http://www.ncbi.nlm.nih.gov/pubmed/16220686)), including the option to cut the features by half each round (instead of one-by-one) if there are many features.

The main function is adapted from <http://www.uccor.edu.ar/paginas/seminarios/Software/SVM_RFE_R_implementation.pdf>.

Also included are tools for wrapping the feature ranking/selection process in an external layer of cross-validation for obtaining unbiased estimates of generalization error/accuracy (See [Ambroise et al., 2002](http://www.ncbi.nlm.nih.gov/pubmed/11983868)).

This entire process can be massively parallelized, and an example implementation is included that uses the [SGE](http://en.wikipedia.org/wiki/Oracle_Grid_Engine) cluster interface.

---
## References

### SVM-RFE
An iterative algorithm that works backward from an initial set of features. At each round it 1) fits a simple linear SVM, 2) ranks the features based on their weights in the SVM solution, and 3) eliminates the feature with the lowest weight.

```
@article{guyon_gene_2002,
    title = {Gene selection for cancer classification using support vector machines},
    volume = {46},
    journal = {Machine learning},
    author = {Guyon, Isabelle and Weston, Jason and Barnhill, Stephen and Vapnik, Vladimir},
    year = {2002},
    note = {Original {SVM-RFE} paper.},
    pages = {389-422}
}
```

### *Multiple* SVM-RFE

Extends this idea by using resampling techniques at each iteration to stabilize the feature rankings. Here we use cross validation. The mSVM-RFE paper is:

```
@article{duan_multiple_2005,
    title = {Multiple SVM-RFE for gene selection in cancer classification with expression data.},
    volume = {4},
    number = {3},
    journal = {IEEE transactions on nanobioscience},
    author = {Duan, Kai-Bo and Rajapakse, Jagath C and Wang, Haiying and Azuaje, Francisco},
    month = sep,
    year = {2005},
    note = {Original Multiple SVM-RFE paper.},
    pages = {228-234}
}
```

---
## Install

1. Download the SVM-RFE package from GitHub.

1. Make sure a recent version of R is available on your computer. http://www.r-project.org

1. Launch R from inside the SVM-RFE directory you just downloaded.

1. Make sure the `e1071` package is installed. This package contains the basic functions for fitting SVMs.

    `install.packages('e1071')`

1. Load the `e1071` package.
    
    `library(e1071)`

1. Source the `msvmRFE.R` document. This makes all the functions available for your R session.
    
    `source('msvmRFE.R')`

---
## Basic usage

For convenience, these commands are compiled in [`demo.R`](https://github.com/johncolby/SVM-RFE/blob/master/demo/demo.R) in the `demo` directory.

### Input

Input data should be formatted as a data frame, with one row for each observation, and one column for each feature. The first column should contain the true class labels. Factor features (ex: gender, site, etc.) should be coded as multiple numeric "dummy" features (see `?model.matrix` for how to automatically generate these in R). 

An example dataset (derived from a bit of the [ADHD200 sample](http://fcon_1000.projects.nitrc.org/indi/adhd200/) is included in the `demo` directory:

```
> set.seed(12345)
> library(e1071)
Loading required package: class
> source('msvmRFE.R')
> load('demo/input.Rdata')
> dim(input)
[1] 208 547
> input[1:5,1:5]
         DX2 Full4.IQ   Age bankssts_SurfArea caudalanteriorcingulate_SurfArea
1000804   TD      109  7.29              1041                              647
1023964 ADHD      123  8.29              1093                              563
1057962 ADHD      129  8.78              1502                              738
1099481 ADHD      116  8.04               826                              475
1127915   TD      124 12.44              1185                              846
```

There are 208 observations (individual subjects), and 546 features. The first column (`DX2`) contains the class labels (`TD`: typically-developing control, `ADHD`: attention deficit hyperactivity disorder).

### svmRFE()

To perform the feature ranking, use the `svmRFE` function:

```
> svmRFE(input, k=10, halve.above=100)
Scaling data...Done!
Features halved from 546 to 273 
Features halved from 273 to 137 
Features halved from 137 to 69 
  |=========================================================================================================| 100%
  [1] 337 319 173 175 503 287 304 492   1  47 286 214 191 402 313 177 187 407 255 468 256 259 391 133 128 113  33
 [28]  68  92 108 497  61 460 383 388 205 465 512 345 107 412 459  60 266 409 504 473 262 263 527  50 520 100 521
 [55] 105 362 172 169 129 532 303  66  85 389  72 370 491 368 365 126 416 400 498 289 218 385  65 131 203 308 542
 [82] 505 396  27 101 485 487 399  22 405 466 212 511 537 359 253  43 111 525 367 392 422 376 186 194 544  90 103
[109] 384  95 149 502 227  89 516 435 463 244  21 395 200 117 298 354 406 248 325 432 144  40 540 479 219 127 420
[136] 440 349 411 423 415 496 217 120 282 452  91 236 301  70 265 478 272 453 535 207 281  58  59 464 543  17  26
[163]  45 397 509 261 331 158 434 297 130  36 356  63 441  41 196 140 403 446 176 472 180 208 268 481  93 123  73
[190] 184 390 104 284 163 450 524 216  23 141  44 486 211 183 518  83 382 166  94 488 311  35 285 546 233 197 273
[217]  31  51 401  80 269  16 274 476  49 160 378 242 121 179  42 136 526 381 161 309 195  24 539 506 327 213 292
[244] 201 188 252 246  34 523 230 146 112 198 353 404 361 226 447 344 501 193 245  52 414 323  48 314 290  46  84
[271] 425 461 342 324 264 182 508 418 449 369 538   7 522 343 189 394 364 455  64 484  86 426 156 437 330 347  11
[298] 408 529 377 295 267 350  28 346 499 371 280 454  74 442  77 185 514 430 142 239 545 489 190 351 145 373 109
[325] 334 102 222 490 375 360 115   8 475 283 114 124 206   2 174 181 277 355  62 240 443 275 387 482 419 519 322
[352]  38   5 451 232 515 224 243 159  56 291 237 321  25  57 316  29 199 221 148 134  30 428 152 513  82 106 341
[379]  67 299 235 448 494  79 234 125 192 410 528 150 366 143 294  76 312 326 357 132 210   6 279  39 424  32 306
[406] 257  13 456  96 170  54 458 530 138 164 474 271 467 433 135 315  87 336 171  10 220 335 534 500 329 457 250
[433] 270 223 483 541 302 168 119 167  37 153   4  88 228 431   9 358 380 469 386 305 155 147 340  19 470 202 333
[460] 328 278 110 293 251  98 151 241 332 317 247  53 288 118 249 495 471 348 531  97 533 137 260 157 231 444  15
[487]  12 238  71 480 258 229 436 307  78 204 154 225 374 427 417 122 116 300 178  81  55 445 310  99  75 536 517
[514] 493 318 338  14 209 462 438 276 507 379 339 139 352 421 254 439  20 413 320 510  18 363 429 372 215 165 477
[541] 162  69 296 398   3 393
```

Here we've indicated that we want `k=10` for the k-fold cross validation as the "multiple" part of mSVM-RFE. To use standard SVM-RFE, you can use `k=1`. Also notice the `halve.above` parameter. This allows you to cut the features in half each round, instead of one by one. This is very useful for data sets with many features. Here we've set `halve.above=100`, so the features will be cut in half each round until there are fewer than 100 remaining. The output is a vector of feature indices, now ordered from most to least "useful".

Note that because of the randomness of the CV draws, features with close ranking scores, and the possible inclusion of useless features, these rankings can change some run to run. However, your output for this demo should be identical because we've all reset the random seed to the same value.

## Estimating generalization error

When exploring machine learning options, it is often useful to estimate generalization error and use this as a benchmark. However, it is important to remember that the **feature selection step must be repeated from scratch on each training set** in whatever cross validation or similar resampling scheme chosen. When feature selection is performed on a data set with many features, it will pick some truly useful features that will generalize, but it will also likely pick some useless features that, by mere chance, happened to align closely with the class labels of the training set. While including these features will give (spuriously) good performance if the error is estimated from this training set itself, the estimated performance will decrease to its true value when the classifier is applied to a true test set where these features are useless. Guyon et al. actually made this mistake in the example demos in their original SVM-RFE paper. This issue is outlined very nicely in:

```
@article{ambroise_selection_2002,
    title = {Selection bias in gene extraction on the basis of microarray gene-expression data.},
    volume = {99},
    number = {10},
    journal = {Proceedings of the National Academy of Sciences of the United States of America},
    author = {Ambroise, Christophe and McLachlan, Geoffrey J},
    month = may,
    year = {2002},
    pages = {6562–6566}
}
```

### Set up folds

Basically, the way to go is to wrap the *entire* feature selection and generalization error estimation process in a top-level loop of external cross validation. For 10-fold CV, we begin by defining which observations are in which folds.

```
> nfold = 10
> nrows = nrow(input)
> folds = rep(1:nfold, len=nrows)[sample(nrows)]
> folds
  [1]  5  8  6 10  9  4  4  1  6  3  2 10  2  6  8  4 10 10  5  8  8  2  6 10  4  5  5  2  2  4  5  9  7  6  2  7
 [37]  4  3  9  2  7  7  4  9  7  1  9  8  5  8  4 10  7  1  1  7  1  6  1  8  1  7  6  9  5  2  5  9  3  3  9  6
 [73]  6  9  2  3  2 10  2  1  6  3  9  6  3  7  5  9  1  8  2  3  7  5  3  6  4  2 10  9 10  8 10  4  3  9  1  4
[109]  7  8  7  1  8  3  8  4  5  7  3  9  8  7  5  3  3  3  5  5  7  2  2  1  9  6  3  1  9  2  4  8  2  5  5  8
[145]  5  9  1  7  4  3  6 10  1 10  8  3  6  4  2 10  3  5  1  8  9  3  7  1  4  1  5  1  3 10 10  1  6  8  4  6
[181]  1  6 10  8  4  6 10  5  5  2  6 10  9  8  4  7  9  4  7  7 10 10  2  8  4  2  6  7
```

In R, many parallel functions like to operate on lists, so next we'll reformat `folds` into a list, where each list element is a vector containing the test set indices for that fold.

```
> folds = lapply(1:nfold, function(x) which(folds == x))
> folds
[[1]]
 [1]   8  46  54  55  57  59  61  80  89 107 112 132 136 147 153 163 168 170 172 176 181

[[2]]
 [1]  11  13  22  28  29  35  40  66  75  77  79  91  98 130 131 138 141 159 190 203 206

[[3]]
 [1]  10  38  69  70  76  82  85  92  95 105 114 119 124 125 126 135 150 156 161 166 173

[[4]]
 [1]   6   7  16  25  30  37  43  51  97 104 108 116 139 149 158 169 179 185 195 198 205

[[5]]
 [1]   1  19  26  27  31  49  65  67  87  94 117 123 127 128 142 143 145 162 171 188 189

[[6]]
 [1]   3   9  14  23  34  58  63  72  73  81  84  96 134 151 157 177 180 182 186 191 207

[[7]]
 [1]  33  36  41  42  45  53  56  62  86  93 109 111 118 122 129 148 167 196 199 200 208

[[8]]
 [1]   2  15  20  21  48  50  60  90 102 110 113 115 121 140 144 155 164 178 184 194 204

[[9]]
 [1]   5  32  39  44  47  64  68  71  74  83  88 100 106 120 133 137 146 165 193 197

[[10]]
 [1]   4  12  17  18  24  52  78  99 101 103 152 154 160 174 175 183 187 192 201 202
```

### Perform feature ranking on all training sets

Using `lapply`, or one of its generic parallel cousins (ex: `sge.parLapply` from the `Rsge` package), we can now perform the feature ranking for all 10 training sets.

```
> results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
...
> length(results)
[1] 10
> results
[[1]]
[[1]]$feature.ids
  [1] 337 173 262  43 486   1 402 308  59 130 535  66 186 505 255  79 354 103 218 542 219 160 504 512 303 537 236
 [28] 273 107  68 391 498 286 412 269 453 266 473 205  33 497  85 409 415 465  89 405 478 187 448  41 133 345 256
 [55]  72  60 422  94 467 450 126 177 502 362 169 384 466 527 491 468 532 392  47 170 265 418 430 368 147 377 282
 [82] 172 492 359 125 411 144  21 101 397 403  22 487 518 113 117 521  17 355 212 420 297 389 464 529 479  50 244
[109] 325 281 321 314  73  91 264 382 216  95  25  27 175 460 180 444 298 503 370 360 324 540  86 419 257 383 433
[136] 263 207 435 108 511 452 331 128 496 349 105 203 140 304 191 152 388 120  65 437  46 356 109  92 396 508 129
[163]  36 522 463 459 440 350 449 319 367 485 148 366 379  13 153  61  90 446 513 517 184 488 131 434 274 193 543
[190] 124 307 171 407  51 253  26 200 293 164 369 390 313 285 116  49 526 201 311 111 230 395 476 353 333 259 159
[217] 385  80 149  28 481 516 217 161 469 318  15 240  99 344 163 346 365 441 232  62 445 295 141 112 197  54 305
[244]  40 351 229 227 183 268 239  71 429 243 316 296 258 104 267 251 252  16 194 179 416  42 380 228 167  37 233
[271] 292 146 461 310  96 393 524 136 400 165 335  76 299 494 339 198 222 458 248 545  83 114 288 408 520 278 414
[298] 100  70 181 348 334  63 500 238 427 399 247 381 134 176 509 470 495  55 279 371 223 115 250  35 320 424 483
[325]  29 214 302 178  23 442 277 533  58 231 138 376  10 249 330  30 208 546 192 462 342 413 404 209 174 501 426
[352] 347 341 394  12 137  44 142 323 290 210 123 332 234 127 195 154 343 361  88 538 340 539 493 471 447 168  52
[379]  81 401 423 306 507 436 155 525 158 294 204 531 428 280 363 241 544 506 534 272 284 150 226 541 338 484 270
[406] 530  74 489 196 514 237 254  64 455  45 245  53  82 221 329 215 166 202  69 271 162  38 289 315 119 283 225
[433] 378 519   9 472  39 242 438 246 189 118 260   2 206   3 185 425  48  24 490 300 182 528 482  84 432 515 261
[460] 312 235 480 474   5  93 454 132 145 375 151 220  56 121 213  32 139 143 287  20  57 475 211   6 406 443  98
[487] 301 456 374 199 417 410 439 364  78 398 336 102 326  67 451  31 188 157 322  75 276 110 317 275 122  11 457
[514]  18 135 357   4 328  77 327 309   8 352 358 477   7 224 373 510 499 291 386  34  87 523 190  19 421  14  97
[541] 431 106 387 372 536 156

[[1]]$train.data.ids
  [1] "1000804" "1023964" "1057962" "1099481" "1127915" "1187766" "1208795" "1320247" "1359325" "1435954" "1471736"
 [12] "1497055" "1511464" "1517240" "1567356" "1700637" "1737393" "1740607" "1780174" "1854959" "1875084" "1884448"
 [23] "1918630" "1934623" "1992284" "1995121" "2030383" "2054438" "2107638" "2136051" "2230510" "2260910" "2306976"
 [34] "2497695" "2570769" "2682736" "2730704" "2735617" "2741068" "2773205" "2821683" "2854839" "2907383" "2950672"
 [45] "2991307" "2996531" "3011311" "3163200" "3174224" "3235580" "3243657" "3433846" "3457975" "3542588" "3619797"
 [56] "3650634" "3653737" "3679455" "3845761" "3999344" "4060823" "4079254" "4084645" "4095229" "4116166" "4154672"
 [67] "4164316" "4187857" "4562206" "4827048" "5164727" "5971050" "6568351" "8009688" "8415034" "8692452" "8697774"
 [78] "8834383" "9326955" "9578663" "9907452" "10001"   "10004"   "10016"   "10009"   "10010"   "10012"   "10018"  
 [89] "10020"   "10021"   "10022"   "10023"   "10028"   "10030"   "10038"   "10034"   "10031"   "10035"   "10037"  
[100] "10039"   "10042"   "10050"   "10051"   "10052"   "10053"   "10056"   "10057"   "10058"   "10059"   "10054"  
[111] "10048"   "10060"   "10047"   "10062"   "10064"   "10065"   "10066"   "10068"   "10069"   "10119"   "10110"  
[122] "10111"   "10112"   "10019"   "10077"   "10078"   "10079"   "10080"   "10008"   "10014"   "10089"   "10090"  
[133] "10091"   "10045"   "10093"   "10094"   "10095"   "10096"   "10098"   "10099"   "10100"   "10101"   "10102"  
[144] "10108"   "10002"   "10011"   "10113"   "10115"   "10006"   "10105"   "10024"   "10116"   "10118"   "10026"  
[155] "10121"   "10122"   "10082"   "10074"   "10075"   "10076"   "10081"   "10083"   "10084"   "10085"   "10036"  
[166] "10086"   "10027"   "10087"   "10025"   "10029"   "10017"   "10049"   "10070"   "10071"   "10033"   "10072"  
[177] "10073"   "10104"   "10088"   "10040"   "10124"   "10125"   "10007"   "10126"   "10127"   "10128"   "10129"  

[[1]]$test.data.ids
 [1] "1283494" "2983819" "3349205" "3349423" "3441455" "3518345" "3601861" "6206397" "9750701" "10032"   "10044"  
[12] "10109"   "10107"   "10092"   "10097"   "10114"   "10106"   "10117"   "10120"   "10123"   "10103"  
...
```

Each list element in `results` contains the feature rankings for that fold (`feature.ids`), as well as the training set row names used to obtain them (`train.data.ids`). The remaining test set row names are included as well (`test.data.ids`).

### Obtain top features across all folds

If we were going to apply these findings to a final test set somewhere, we would still want the best features across *all* of this training data.

```
> top.features = WriteFeatures(results, input, save=F)
> head(top.features)
                      FeatureName FeatureID AvgRank
1     superiortemporal_MeanCurv.1       337     2.8
2             bankssts_ThickAvg.1       173     5.6
3   posteriorcingulate_ThickStd.1       262     7.9
4 rostralmiddlefrontal_GausCurv.1       402    26.0
5            frontalpole_SurfArea        33    27.2
6         temporalpole_SurfArea.1        68    37.2
```

Ordered by average rank across the 10 folds (`AvgRank`, lower numbers are better), this gives us a list of the feature names (`FeatureName`, i.e. the corresponding column name from `input`), as well as the feature indices (`FeatureID`, i.e. the corresponding column index from `input` minus 1).

### Estimate generalization error using a varying number of top features

Now that we have a ranking of features for each of the 10 training sets, the final step is to estimate the generalization error we can expect if we were to train a final classifier on these features and apply it to a new test set. Here, a radial basis function kernel SVM is tuned on each training set independently. This consists of doing internal 10-fold CV error estimation at each combination of SVM hyperparameters (Cost and Gamma) using grid search. The optimal parameters are then used to train the SVM on the entire training set. Finally, generalization error is determined by predicting the corresponding test set. This is done for each fold in the external 10-fold CV, and all 10 of these generalization error estimates are averaged to give more stability. This process is repeated while varying the number of top features that are used as input, and there will typically be a "sweet spot" where there are not too many nor too few features. Outlined, this process looks like:

```
external 10x CV
  Rank features with mSVM-RFE
  for(nfeat=1 to 500ish)
    Grid search over SVM parameters
      10x CV
        Train SVM
        Obtain generalization error estimate
      Average generalization error estimate across multiple folds
    Choose parameters with best average performance
    Train SVM on full training set
    Obtain generalization error estimate on corresponding external CV test set
Average generalization errors across multiple folds
Choose the optimum number of features
```

To implement it over the top 5 features, we do:

```
> featsweep = lapply(1:5, FeatSweep.wrap, results, input)
> featsweep
[[1]]
[[1]]$svm.list
[[1]]$svm.list[[1]]
   gamma cost     error dispersion
1 0.0625   16 0.3333333         NA

[[1]]$svm.list[[2]]
    gamma cost     error dispersion
1 0.03125    4 0.2857143         NA

[[1]]$svm.list[[3]]
  gamma cost     error dispersion
1 0.125    4 0.6666667         NA
...

[[1]]$error
[1] 0.5290476

... ...
```

Each `featsweep` list element corresponds to using that many of the top features (i.e. `featsweep[1]` is using only the top feature, `featsweep[2]` is using the top 2 features, etc.). Within each, `svm.list` contains the generalization error estimates for each of the 10 folds in the external 10-fold CV. These accuracies are averaged as `error`.

### Plot of generalization error vs. # of features

To show these results visually, we can plot the average generalization error vs. the number of top features used as input. 

For reference, it is useful to show the chance error rate. Typically, this is equal to the "no information" rate we would get if we simply always picked the class label with the greater prevalence in the data set.

```
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info)
dev.off()
```

![msvmrfe_demo](https://user-images.githubusercontent.com/473295/88239232-a9f71a80-cc38-11ea-8898-d43864d1d0ff.png)

## Parallelization

As you can probably see, the main limitation in doing this type of exploration is processing time. For example, in this demonstration, and just considering the number of times we have to fit an SVM, we have:

- **Feature ranking**

    10 external folds x 546 features x 10 msvmRFE folds = 54600 linear SVMs

- **Generalization error estimation**

    10 external folds x 546 features x 169 hyperparamter combos for exhaustive search x 10 folds each = 9227400 RBF kernel SVMs

We have already shortened this some by 1) eliminating more than one feature at a time in the feature ranking step, and 2) only estimating generalization accuracies across a sweep of the top features.

Other ideas to shorten processing time:

- Fewer external CV folds.
- Smaller or coarser grid for parameter sweep in SVM tuning step
- Fewer CV folds in SVM tuning step.

This code is already set up to use `lapply` calls for these 2 mains task, so fortunately, they can be relatively easily parallelized using a variety of R packages that include parallel versions of `lapply`. Examples include:

- `mclapply` in the [`multicore`](http://cran.r-project.org/web/packages/multicore/index.html) package
- `sge.parLapply` in the [`Rsge`](http://cran.r-project.org/web/packages/Rsge/index.html) package

## Example SGE implementation

Included in the `sge` directory is an example custom parallel implementation using the SGE cluster interface. It will need to be modified to fit your own cluster setup, but I included it because looking through the code might give you some ideas.

- For the **feature ranking**, the external CV folds are parallelized.
- For the **generalization error estimating**, all the different runs with varying numbers of input features are parallelized.

### Setup notes

- A version of `Rscript` must be in your unix path. `Rscript` can be found in the `bin` directory of recent versions of R.
- These scripts are executed like shell scripts, so...
  - Make sure they are in your unix path (e.g. something like `PATH=/path/to/SVM-RFE/sge:$PATH` to add the `sge` directory to the path in bash).
  - Make sure they are executable (e.g. `chmod +x sge/*`).
- As written, these scripts use an excellent `qsub` wrapper called `fsl_sub` that is included with the [FSL](http://www.fmrib.ox.ac.uk/fsl/) neuroimaging toolset. If yours is different, you'll need to modify these calls accordingly.
- If `e1071` is installed, but there are error messages that it can't be found, you probably installed it in some user specific location. You can add this area to the default R search path by setting the `R_LIBS` unix environment variable.
- The main `msvmRFE.R` functions file is sourced in several of these scripts. It is setup to work relative to the `demo` directory, as downloaded, but you might want to hard-code your system-specific full path to this file so you don't always have to call it from the same directory.

### Usage

From within a directory with an `input.Rdata` file present, submit the feature ranking jobs from the unix command line.

```
$ msvmRFE_submit.R nfold=10
1484221
1484222
Number of folds: 10

```

You can make sure your jobs are running smoothly with `qstat`.

```
$ watch qstat -u me
Every 2.0s: qstat -u me                                                                Wed Sep 14 10:30:03 2011

job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
1484221 0.50093 jobarray.t me           r     09/14/2011 11:20:44 pod_smp.q@n7261                    1 1
1484221 0.50046 jobarray.t me           r     09/14/2011 11:20:44 pod_smp.q@n6266                    1 2
1484221 0.50031 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n6228                    1 3
1484221 0.50023 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n6231                    1 4
1484221 0.50018 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n7225                    1 5
1484221 0.50015 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n6273                    1 6
1484221 0.50013 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n7223                    1 7
1484221 0.50012 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n6237                    1 8
1484221 0.50010 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n7231                    1 9
1484221 0.50009 jobarray.t me           r     09/14/2011 11:20:45 pod_smp.q@n6233                    1 10
1484222 0.00000 msvmRFE_co me           hqw   09/14/2011 11:19:42                                    1        
```

As expected, the 10 external cross validation folds are running in parallel, and 1 extra job is waiting to combine those results.

You can also check the progress of these jobs.

```
$ tail jobarray.txt.o*
==> jobarray.txt.o1484221.1 <==
Scaling data...Done!
  |===================================                                   |  49%
==> jobarray.txt.o1484221.10 <==
Scaling data...Done!
  |================================                                      |  46%
==> jobarray.txt.o1484221.2 <==
Scaling data...Done!
  |===================================                                   |  50%
==> jobarray.txt.o1484221.3 <==
Scaling data...Done!
  |========================================                              |  57%
==> jobarray.txt.o1484221.4 <==
Scaling data...Done!
  |================================                                      |  45%
==> jobarray.txt.o1484221.5 <==
Scaling data...Done!
  |=============================                                         |  41%
==> jobarray.txt.o1484221.6 <==
Scaling data...Done!
  |======================================                                |  54%
==> jobarray.txt.o1484221.7 <==
Scaling data...Done!
  |===================================                                   |  49%
==> jobarray.txt.o1484221.8 <==
Scaling data...Done!
  |=====================================                                 |  53%
==> jobarray.txt.o1484221.9 <==
Scaling data...Done!
  |======================================                                |  55%
```

Once these jobs have finished, you'll be left with a `results.Rdata` file, containing the results of the feature ranking for each of the external CV folds (like the interactive example, above), as well as a plain text version of the average feature rankings across all folds (`features_ranked.txt`).

Next, submit the generalization error estimation sweep over the top features.

```
$ featsweep_submit.R nfeat=1:5
1484244
1484245
Number of features to sweep: 5
```

Input arguments to Rscript get parsed as normal R expressions, so you can request a more complicated set if desired. For example:

```
$ featsweep_submit.R 'nfeat=c(1:100, seq(102, 500, by=2))'
```

When these jobs finish, you should have a `featsweep.Rdata` file, containing the obtimal tuning parameters and generalization error estimates for each combination of top features (like the interactive example, above), as well as 2 plots of the estimated generalization error vs. # of top features. 