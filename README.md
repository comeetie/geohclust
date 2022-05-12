
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geohclust

<!-- badges: start -->
<!-- badges: end -->

geohclust offers two functions `?geohclust_poly` and `?geohclust_graph`
that enable the clustering of spatial data such as polygons with a
hclust type approach but taking advantages of contiguity constraints.
The contiguity naturally create a sparsely connected graph that can be
leveraged to speed-up the calculations and deal with more than 30000
polygons in seconds.

## Installation

You can install the development version of geohclust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("comeetie/geohclust")
```

## Example

This is a basic example, we first prepare some spatial polygons data,
here the results at the municipality level in one french department for
the :

``` r
library(geohclust)
library(dplyr)
library(sf)
data("modesshare")
deplist = c(37)
dep = modesshare |> 
  filter(DEP %in% deplist) 
```

Do the clustering and use the classical function from `?hclust`
(`?plot.hclust` and `?cutree`):

``` r
hc=geohclust_poly(dep)
#> Warning: Some features were not numeric and have been removed from the
#> clustering.
#> Merge :174--212
#> Merge :349--364
#> Merge :105--193
#> Merge :90--243
#> Merge :44--238
#> Merge :94--358
#> Merge :68--117
#> Merge :35--102
#> Merge :56--253
#> Merge :111--203
#> Merge :307--322
#> Merge :285--383
#> Merge :339--379
#> Merge :195--320
#> Merge :65--225
#> Merge :53--385
#> Merge :72--377
#> Merge :244--318
#> Merge :141--211
#> Merge :310--378
#> Merge :146--370
#> Merge :58--209
#> Merge :129--297
#> Merge :178--401
#> Merge :153--283
#> Merge :239--382
#> Merge :165--404
#> Merge :222--348
#> Merge :162--223
#> Merge :175--381
#> Merge :25--121
#> Merge :218--252
#> Merge :39--132
#> Merge :130--308
#> Merge :182--247
#> Merge :16--354
#> Merge :98--326
#> Merge :47--269
#> Merge :1--356
#> Merge :22--294
#> Merge :271--374
#> Merge :258--415
#> Merge :76--337
#> Merge :54--387
#> Merge :389--408
#> Merge :55--184
#> Merge :287--363
#> Merge :233--380
#> Merge :273--426
#> Merge :275--392
#> Merge :96--399
#> Merge :26--344
#> Merge :113--126
#> Merge :120--321
#> Merge :143--416
#> Merge :136--312
#> Merge :125--352
#> Merge :304--420
#> Merge :147--202
#> Merge :234--259
#> Merge :31--295
#> Merge :319--439
#> Merge :342--367
#> Merge :397--419
#> Merge :421--431
#> Merge :27--95
#> Merge :208--293
#> Merge :205--309
#> Merge :40--290
#> Merge :194--447
#> Merge :436--448
#> Merge :46--280
#> Merge :296--317
#> Merge :314--406
#> Merge :272--279
#> Merge :103--291
#> Merge :69--407
#> Merge :219--331
#> Merge :232--353
#> Merge :108--405
#> Merge :67--458
#> Merge :78--390
#> Merge :45--398
#> Merge :231--369
#> Merge :29--114
#> Merge :384--393
#> Merge :180--422
#> Merge :281--412
#> Merge :190--423
#> Merge :84--116
#> Merge :357--388
#> Merge :118--323
#> Merge :257--445
#> Merge :52--303
#> Merge :64--207
#> Merge :395--403
#> Merge :254--361
#> Merge :77--133
#> Merge :140--454
#> Merge :86--386
#> Merge :333--424
#> Merge :164--432
#> Merge :413--446
#> Merge :48--189
#> Merge :249--262
#> Merge :93--151
#> Merge :115--334
#> Merge :226--484
#> Merge :366--418
#> Merge :8--325
#> Merge :80--210
#> Merge :75--127
#> Merge :428--460
#> Merge :109--491
#> Merge :270--429
#> Merge :455--459
#> Merge :59--181
#> Merge :214--375
#> Merge :74--316
#> Merge :206--242
#> Merge :61--498
#> Merge :346--461
#> Merge :89--400
#> Merge :154--311
#> Merge :449--464
#> Merge :13--503
#> Merge :236--504
#> Merge :79--476
#> Merge :159--302
#> Merge :81--168
#> Merge :329--450
#> Merge :467--492
#> Merge :18--505
#> Merge :188--453
#> Merge :9--480
#> Merge :145--191
#> Merge :32--360
#> Merge :434--451
#> Merge :442--508
#> Merge :148--517
#> Merge :173--266
#> Merge :306--444
#> Merge :215--402
#> Merge :245--263
#> Merge :217--522
#> Merge :131--409
#> Merge :494--524
#> Merge :171--525
#> Merge :425--485
#> Merge :150--435
#> Merge :350--469
#> Merge :73--362
#> Merge :228--376
#> Merge :251--282
#> Merge :457--468
#> Merge :134--166
#> Merge :37--471
#> Merge :2--466
#> Merge :7--85
#> Merge :355--373
#> Merge :220--313
#> Merge :301--470
#> Merge :196--204
#> Merge :42--541
#> Merge :24--443
#> Merge :138--490
#> Merge :28--237
#> Merge :520--545
#> Merge :107--330
#> Merge :149--535
#> Merge :41--438
#> Merge :427--437
#> Merge :221--343
#> Merge :529--551
#> Merge :0--414
#> Merge :324--475
#> Merge :265--359
#> Merge :152--365
#> Merge :169--548
#> Merge :298--335
#> Merge :248--558
#> Merge :51--559
#> Merge :34--560
#> Merge :160--440
#> Merge :163--465
#> Merge :106--276
#> Merge :101--532
#> Merge :5--60
#> Merge :158--502
#> Merge :63--477
#> Merge :300--463
#> Merge :198--229
#> Merge :338--570
#> Merge :192--345
#> Merge :88--521
#> Merge :507--568
#> Merge :472--537
#> Merge :157--478
#> Merge :235--274
#> Merge :14--574
#> Merge :3--119
#> Merge :268--527
#> Merge :340--550
#> Merge :433--581
#> Merge :10--144
#> Merge :128--526
#> Merge :410--584
#> Merge :394--585
#> Merge :347--518
#> Merge :486--513
#> Merge :110--411
#> Merge :336--589
#> Merge :497--590
#> Merge :288--292
#> Merge :123--289
#> Merge :327--487
#> Merge :30--177
#> Merge :17--499
#> Merge :230--260
#> Merge :36--516
#> Merge :511--598
#> Merge :33--62
#> Merge :66--536
#> Merge :372--599
#> Merge :91--602
#> Merge :139--571
#> Merge :104--368
#> Merge :538--593
#> Merge :315--549
#> Merge :176--483
#> Merge :456--605
#> Merge :479--609
#> Merge :569--573
#> Merge :70--156
#> Merge :299--523
#> Merge :200--277
#> Merge :50--531
#> Merge :199--603
#> Merge :20--563
#> Merge :87--305
#> Merge :227--351
#> Merge :38--586
#> Merge :240--543
#> Merge :510--621
#> Merge :501--576
#> Merge :250--542
#> Merge :43--441
#> Merge :554--625
#> Merge :417--588
#> Merge :530--627
#> Merge :371--506
#> Merge :328--623
#> Merge :186--630
#> Merge :12--607
#> Merge :582--622
#> Merge :452--633
#> Merge :587--634
#> Merge :241--628
#> Merge :509--636
#> Merge :528--637
#> Merge :246--579
#> Merge :515--639
#> Merge :170--606
#> Merge :391--620
#> Merge :617--642
#> Merge :112--619
#> Merge :97--591
#> Merge :557--567
#> Merge :601--646
#> Merge :519--641
#> Merge :135--648
#> Merge :92--496
#> Merge :286--649
#> Merge :99--644
#> Merge :597--632
#> Merge :604--647
#> Merge :652--654
#> Merge :572--624
#> Merge :6--481
#> Merge :474--645
#> Merge :19--655
#> Merge :615--659
#> Merge :534--613
#> Merge :187--495
#> Merge :657--662
#> Merge :430--663
#> Merge :512--596
#> Merge :264--665
#> Merge :185--556
#> Merge :578--667
#> Merge :267--668
#> Merge :82--631
#> Merge :561--594
#> Merge :278--653
#> Merge :500--658
#> Merge :616--643
#> Merge :255--610
#> Merge :552--670
#> Merge :493--676
#> Merge :566--674
#> Merge :201--626
#> Merge :675--679
#> Merge :172--660
#> Merge :540--592
#> Merge :547--682
#> Merge :473--564
#> Merge :462--664
#> Merge :583--661
#> Merge :553--678
#> Merge :635--687
#> Merge :514--688
#> Merge :256--689
#> Merge :539--690
#> Merge :124--666
#> Merge :488--575
#> Merge :546--677
#> Merge :562--694
#> Merge :396--681
#> Merge :167--577
#> Merge :555--669
#> Merge :638--698
#> Merge :57--691
#> Merge :629--693
#> Merge :533--684
#> Merge :672--701
#> Merge :71--261
#> Merge :15--600
#> Merge :611--700
#> Merge :224--673
#> Merge :183--699
#> Merge :686--703
#> Merge :595--697
#> Merge :482--702
#> Merge :544--711
#> Merge :179--706
#> Merge :197--341
#> Merge :580--650
#> Merge :565--695
#> Merge :608--696
#> Merge :708--717
#> Merge :618--683
#> Merge :216--692
#> Merge :23--720
#> Merge :489--713
#> Merge :716--722
#> Merge :704--710
#> Merge :83--718
#> Merge :671--723
#> Merge :142--725
#> Merge :213--727
#> Merge :332--612
#> Merge :21--680
#> Merge :640--709
#> Merge :651--724
#> Merge :4--728
#> Merge :712--719
#> Merge :284--734
#> Merge :49--732
#> Merge :685--731
#> Merge :656--737
#> Merge :100--730
#> Merge :161--726
#> Merge :705--736
#> Merge :137--735
#> Merge :707--721
#> Merge :740--743
#> Merge :715--733
#> Merge :614--745
#> Merge :714--741
#> Merge :742--747
#> Merge :122--744
#> Merge :738--746
#> Merge :739--749
#> Merge :729--750
#> Merge :155--748
#> Merge :11--751
#> Merge :752--754
#> Merge :753--755
plot(hc)
cutree(hc,k=30) |> head(20)
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  2  2  3  4  1  5  6  6  2  6  7  6  1  2  8  1  9  1  2
```

<img src="man/figures/README-clustering-1.png" width="100%" />

You may also use the `?geocutree` function which build directly a
spatial data.frame with the clustering results:

``` r
plot(geocutree(hc,k=30))
```

<img src="man/figures/README-plot-1.png" width="100%" />
