
# Accessing strain survey data using MPD

In this example we will programmatically retrieve a measure record from MPD. The measure we are working with here reflects total locomotor activity 3h post cocaine injection over a test period of 4 hours. This measure has an MPD measure ID of 40708 and is named activ_cocaine2_3. A total of 15 strains were assayed across both sexs with an average n of 8 mice per sex per strain. 

We will first retrieve individual animal values for this measure.


```python
import requests
import json
from pandas.io.json import json_normalize

resp = requests.get("https://phenome.jax.org/api/pheno/animalvals/40708") 
data = json.loads(resp.text)
animaldata = json_normalize(data['animaldata'])
print (animaldata)

#with open("queryMPD_SingleMeasure.tsv", 'w') as fl:    
#    fl.writelines(resp.text)
```

        animal_id  measnum   projsym sex stocknum         strain  strainid  \
    0         639    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    1         638    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    2         637    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    3         636    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    4         643    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    5         642    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    6         641    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    7         640    40708  Thomsen1   f   002448    129S1/SvImJ         3   
    8         631    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    9         633    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    10        634    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    11        635    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    12        628    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    13        629    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    14        630    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    15        632    40708  Thomsen1   m   002448    129S1/SvImJ         3   
    16        753    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    17        752    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    18        751    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    19        750    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    20        749    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    21        748    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    22        747    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    23        754    40708  Thomsen1   f      809  129S6/SvEvTac        53   
    24        745    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    25        739    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    26        740    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    27        741    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    28        742    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    29        743    40708  Thomsen1   m      809  129S6/SvEvTac        53   
    ..        ...      ...       ...  ..      ...            ...       ...   
    202       755    40708  Thomsen1   m   000686          SJL/J        17   
    203       655    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    204       657    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    205       651    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    206       652    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    207       653    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    208       654    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    209       658    40708  Thomsen1   f   001146      SPRET/EiJ        39   
    210       644    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    211       646    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    212       647    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    213       650    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    214       649    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    215       648    40708  Thomsen1   m   001146      SPRET/EiJ        39   
    216       486    40708  Thomsen1   f      856         Tac:SW        54   
    217       489    40708  Thomsen1   f      856         Tac:SW        54   
    218       490    40708  Thomsen1   f      856         Tac:SW        54   
    219       491    40708  Thomsen1   f      856         Tac:SW        54   
    220       492    40708  Thomsen1   f      856         Tac:SW        54   
    221       487    40708  Thomsen1   f      856         Tac:SW        54   
    222       488    40708  Thomsen1   f      856         Tac:SW        54   
    223       485    40708  Thomsen1   f      856         Tac:SW        54   
    224       484    40708  Thomsen1   m      856         Tac:SW        54   
    225       483    40708  Thomsen1   m      856         Tac:SW        54   
    226       482    40708  Thomsen1   m      856         Tac:SW        54   
    227       481    40708  Thomsen1   m      856         Tac:SW        54   
    228       480    40708  Thomsen1   m      856         Tac:SW        54   
    229       479    40708  Thomsen1   m      856         Tac:SW        54   
    230       478    40708  Thomsen1   m      856         Tac:SW        54   
    231       477    40708  Thomsen1   m      856         Tac:SW        54   
    
          value           varname  zscore  
    0     497.0  activ_cocaine2_3   -0.08  
    1     631.0  activ_cocaine2_3    0.87  
    2     609.0  activ_cocaine2_3    0.72  
    3     264.0  activ_cocaine2_3   -1.74  
    4     699.0  activ_cocaine2_3    1.36  
    5     488.0  activ_cocaine2_3   -0.15  
    6     382.0  activ_cocaine2_3   -0.90  
    7     499.0  activ_cocaine2_3   -0.07  
    8     452.0  activ_cocaine2_3    0.05  
    9     371.0  activ_cocaine2_3   -0.50  
    10    236.0  activ_cocaine2_3   -1.42  
    11    239.0  activ_cocaine2_3   -1.39  
    12    573.0  activ_cocaine2_3    0.87  
    13    555.0  activ_cocaine2_3    0.75  
    14    528.0  activ_cocaine2_3    0.56  
    15    604.0  activ_cocaine2_3    1.08  
    16     31.0  activ_cocaine2_3   -2.14  
    17    488.0  activ_cocaine2_3    0.44  
    18    374.0  activ_cocaine2_3   -0.20  
    19    502.0  activ_cocaine2_3    0.52  
    20    480.0  activ_cocaine2_3    0.39  
    21    487.0  activ_cocaine2_3    0.43  
    22    609.0  activ_cocaine2_3    1.12  
    23    309.0  activ_cocaine2_3   -0.57  
    24    126.0  activ_cocaine2_3   -1.36  
    25    574.0  activ_cocaine2_3    1.33  
    26    452.0  activ_cocaine2_3    0.60  
    27    561.0  activ_cocaine2_3    1.25  
    28    190.0  activ_cocaine2_3   -0.97  
    29    322.0  activ_cocaine2_3   -0.18  
    ..      ...               ...     ...  
    202   534.0  activ_cocaine2_3   -0.94  
    203  1065.0  activ_cocaine2_3    1.78  
    204   597.0  activ_cocaine2_3   -0.37  
    205   474.0  activ_cocaine2_3   -0.94  
    206   488.0  activ_cocaine2_3   -0.87  
    207   875.0  activ_cocaine2_3    0.91  
    208   575.0  activ_cocaine2_3   -0.47  
    209   673.0  activ_cocaine2_3   -0.02  
    210   821.0  activ_cocaine2_3    1.07  
    211   807.0  activ_cocaine2_3    1.02  
    212   652.0  activ_cocaine2_3    0.46  
    213   382.0  activ_cocaine2_3   -0.50  
    214   336.0  activ_cocaine2_3   -0.66  
    215   134.0  activ_cocaine2_3   -1.39  
    216  2571.0  activ_cocaine2_3   -0.53  
    217  2868.0  activ_cocaine2_3   -0.26  
    218  2217.0  activ_cocaine2_3   -0.87  
    219  2932.0  activ_cocaine2_3   -0.20  
    220  2203.0  activ_cocaine2_3   -0.88  
    221  2983.0  activ_cocaine2_3   -0.15  
    222  3950.0  activ_cocaine2_3    0.76  
    223  5407.0  activ_cocaine2_3    2.12  
    224   776.0  activ_cocaine2_3   -0.61  
    225  1071.0  activ_cocaine2_3   -0.47  
    226  1578.0  activ_cocaine2_3   -0.24  
    227  2636.0  activ_cocaine2_3    0.26  
    228   943.0  activ_cocaine2_3   -0.53  
    229  2196.0  activ_cocaine2_3    0.05  
    230   421.0  activ_cocaine2_3   -0.78  
    231  7040.0  activ_cocaine2_3    2.32  
    
    [232 rows x 10 columns]


# Assess strain x sex effects

Now let us fit a linear model that assesses the main effects of strain and sex in addition to strain by sex effects



```python
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.graphics.factorplots import interaction_plot

formula = 'value ~ C(sex) + C(strain) + C(strain):C(sex)'
model = ols(formula, animaldata).fit()
aov_table = anova_lm(model)

print(aov_table)

```

                         df        sum_sq       mean_sq          F        PR(>F)
    C(sex)              1.0  2.931693e+07  2.931693e+07   8.435830  4.088519e-03
    C(strain)          14.0  1.172687e+09  8.376333e+07  24.102563  9.216382e-36
    C(strain):C(sex)   14.0  1.231190e+08  8.794214e+06   2.530500  2.365076e-03
    Residual          202.0  7.020080e+08  3.475287e+06        NaN           NaN


# Estimating heritability

Quite simply, heritability is defined as the ratio of genotypic and phenotypic variance. Therefore to estimate heritability we will need to obtain the between subject variance (genetic variance) along with the within subject variance (error variance) and interaction variance contained in the strain:sex term.


```python
msb = aov_table.mean_sq[1] #between subject variance
mse = aov_table.mean_sq[3] #within subject variance
msi = aov_table.mean_sq[2] #interaction variance
s = 2
n = 8

h2 = round((msb-msi)/(s*(msb+(n-1)*mse)),2)

print("heritability for this measure is : ", h2)

```

    heritability for this measure is :  0.35

