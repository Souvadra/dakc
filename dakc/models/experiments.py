# this notebook defines all the input parameters inside each synthetic datasets

inputlist = ["s20", "s21", "s22", "s23", "s24", "s25", "s26", "s27", "s28", "s29", "s30"] #, "s31", "s32"]

num_reads = {}
num_reads["s20"] = 349500
num_reads["s21"] = 699050
num_reads["s22"] = 1398100
num_reads["s23"] = 2796200
num_reads["s24"] = 5592400
num_reads["s25"] = 11184800
num_reads["s26"] = 22369600
num_reads["s27"] = 44739200
num_reads["s28"] = 89478450
num_reads["s29"] = 178956950
num_reads["s30"] = 357913900
# num_reads["s31"] = 715827850
# num_reads["s32"] = 1431655750

total_time = {}
total_time["s20"] = 0.025
total_time["s21"] = 0.04
total_time["s22"] = 0.073
total_time["s23"] = 0.134
total_time["s24"] = 0.25
total_time["s25"] = 0.48
total_time["s26"] = 0.99
total_time["s27"] = 2.00
total_time["s28"] = 4.11
total_time["s29"] = 8.61
total_time["s30"] = 20.04

p1_time = {}
p1_time["s20"] = 0.02
p1_time["s21"] = 0.03
p1_time["s22"] = 0.043
p1_time["s23"] = 0.089
p1_time["s24"] = 0.145
p1_time["s25"] = 0.28
p1_time["s26"] = 0.54 
p1_time["s27"] = 1.12
p1_time["s28"] = 2.13
p1_time["s29"] = 4.6
p1_time["s30"] = 9.3

p2_time = {}
for input in inputlist:
  p2_time[input] = total_time[input] - p1_time[input]

sort_cache = {}
sort_cache["s20"] = 1213856
sort_cache["s21"] = 4417528
sort_cache["s22"] = 11263584
sort_cache["s23"] = 24833520
sort_cache["s24"] = 52136376
sort_cache["s25"] = 107445312
sort_cache["s26"] = 251671672
sort_cache["s27"] = 597555616
sort_cache["s28"] = 1306281240
sort_cache["s29"] = 2729037576
sort_cache["s30"] = 5581262928

scache_err = {}
scache_err["s20"] = 9440.59
scache_err["s21"] = 22322.86
scache_err["s22"] = 8252.54
scache_err["s23"] = 21170.22
scache_err["s24"] = 13240.06
scache_err["s25"] = 7152.36
scache_err["s26"] = 260638.67
scache_err["s27"] = 839887.67
scache_err["s28"] = 641482.51
scache_err["s29"] = 678586.39
scache_err["s30"] = 1251367.46

parse_cache = {}
parse_cache["s20"] = 316079.30        # 118400
parse_cache["s21"] = 581177.62        # 226896
parse_cache["s22"] = 1072155.58       # 436536
parse_cache["s23"] = 5307680.042      # 856512
parse_cache["s24"] = 13367847.79      # 1689096
parse_cache["s25"] = 29334681.50      # 3350960
parse_cache["s26"] = 61609726.04      # 6680056
parse_cache["s27"] = 125434034.12     # 13324344
parse_cache["s28"] = 253610567.71     # 26624208
parse_cache["s29"] = 510404030.17     # 53204416
parse_cache["s30"] = 1024381690.21    # 106327016

pcache_err = {}
pcache_err["s20"] = 3621.44
pcache_err["s21"] = 1197.45
pcache_err["s22"] = 22256.67
pcache_err["s23"] = 37250.30
pcache_err["s24"] = 30996.52
pcache_err["s25"] = 162381.61
pcache_err["s26"] = 211257.62
pcache_err["s27"] = 79898.34
pcache_err["s28"] = 201092.323
pcache_err["s29"] = 913257.399
pcache_err["s30"] = 1399211.54
