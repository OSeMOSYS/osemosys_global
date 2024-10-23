"""Constants for the transmission module"""

SET_DTYPES = {
    "DAILYTIMEBRACKET": int,
    "EMISSION":str,
    "FUEL":str,
    "MODE_OF_OPERATION":int,
    "REGION":str,
    "SEASON":str,
    "STORAGE":str,
    "TECHNOLOGY":str,
    "TIMESLICE":str,
    "YEAR":int,
}

"""Set the year at which existing transmission capacity is 
assumed to be retired."""
RETIREMENT_YEAR_TRANSMISSION = 2060

"""Set the build year for planned transmission capacity entries
from the GTD dataset that do not have a expected commissioning 
year attached."""
PLANNED_BUILD_YEAR_TRANSMISSION = 2030

"""Set the 'From' OG region to which the residual capacity
based on the GTD dataset entry should be allocated. E.g. USAMI 
in the GTD dataset consists of multiple OG regions (e.g. USAME,
USASW, USASA) but the residual capacity has to be allocated
to one individual region."""
CUSTOM_TRN_BA_DICT_FROM = {
    'TRNUSAMIUSAPJ' : 'USASW',
    'TRNUSAMIUSASO' : 'USASA',
    'TRNUSAMIUSASP' : 'USASA',
    'TRNUSAMIUSASW' : 'USASW',
    'TRNUSAMIUSATV' : 'USASA',
    'TRNUSAMIUSAWA' : 'USAMW',
    'TRNUSAPJUSATV' : 'USARW',
    'TRNUSASWUSAWA' : 'USAMW',
    'TRNUSASWUSAWM' : 'USASN',
    'TRNCHNNMCHNSD' : 'CHNEM',
    'TRNCHNNMCHNSI' : 'CHNWM',
    'TRNCHNNMCHNSX' : 'CHNWM',
    'TRNCHNNMCHNTJ' : 'CHNEM',
    'TRNCHNNMMNGXX' : 'CHNWM',
    'TRNCHNNMRUSSI' : 'CHNEM',
}

"""Set the 'To' OG region to which the residual capacity
based on the GTD dataset entry should be allocated. E.g. USAMI 
in the GTD dataset consists of multiple OG regions (e.g. USAME,
USASW, USASA) but the residual capacity has to be allocated
to one individual region."""
CUSTOM_TRN_BA_DICT_TO = {
    'TRNCANMBUSAMI' : 'USAMW',
    'TRNCANMBUSASW' : 'USAMW',
    'TRNCANONUSAMI' : 'USARM',
    'TRNCANONUSAPJ' : 'USARE',
    'TRNCANSKUSAMI' : 'USAMW',
    'TRNCANSKUSASW' : 'USAMW',
    'TRNUSAAIUSAMI' : 'USASW',
    'TRNUSAAIUSASW' : 'USASN',
    'TRNUSACEUSAPJ' : 'USARE',
    'TRNUSACWUSAPJ' : 'USARW',
    'TRNUSADUUSAPJ' : 'USARW',
    'TRNUSAEPUSASW' : 'USASS',
    'TRNUSAERUSAMI' : 'USASA',
    'TRNUSAERUSASW' : 'USASS',
    'TRNUSAGLUSAMI' : 'USASW',
    'TRNUSALGUSAMI' : 'USASW',
    'TRNUSALGUSAPJ' : 'USARW',
    'TRNUSAMIUSAPJ' : 'USARW',
    'TRNUSAMIUSASW' : 'USASN',
    'TRNUSANYUSAPJ' : 'USARE',
    'TRNUSAPNUSASW' : 'USASS',
    'TRNUSAPOUSASW' : 'USASN',
    'TRNUSASPUSASW' : 'USASS',
    'TRNCHNBECHNNM' : 'CHNWM',
    'TRNCHNGACHNNM' : 'CHNWM',
    'TRNCHNHBCHNNM' : 'CHNEM',
    'TRNCHNHJCHNNM' : 'CHNEM',
    'TRNCHNJICHNNM' : 'CHNEM',
    'TRNCHNJXCHNNM' : 'CHNWM',
    'TRNCHNLICHNNM' : 'CHNEM',
    'TRNCHNNICHNNM' : 'CHNWM',
}

"""Set the OG transmission technologies that are missing in the 
GTD dataset."""
CUSTOM_TRN_BA_MISSING = [
    'TRNCANONUSAMW',
    'TRNCHNBECHNEM',
    'TRNCHNEMCHNWM',
    'TRNCHNEMMNGXX',
    'TRNCHNHBCHNWM',
    'TRNUSAMEUSAMW',
    'TRNUSAMEUSARW',
    'TRNUSAMWUSARA',
    'TRNUSAMWUSARW',
    'TRNUSAMWUSASN',
    'TRNUSAMWUSASW',
    'TRNUSARAUSASS',
    'TRNUSAREUSARW',
    'TRNUSAREUSASV',
    'TRNUSASAUSASN',
    'TRNUSASAUSASW',
    'TRNUSASNUSASS',    
    ]

"""Set transmission pathways that should be based on HVDC_subsea."""
SUBSEA_LINES = [
    'TRNAREXXINDWE',
    'TRNAREXXIRNXX',
    'TRNAUSQLPNGXX',
    'TRNAUSTAAUSVI',
    'TRNAUSWAIDNNU',
    'TRNAUSWATLSXX',
    'TRNBELXXGBRXX',
    'TRNCANARCANNL',
    'TRNCANNLGRLXX',
    'TRNCANONUSARE',
    'TRNCHNFUTWNXX',
    'TRNCHNSDKORXX',
    'TRNCHNSDPRKXX',
    'TRNCYPXXEGYXX',
    'TRNCYPXXGRCXX',
    'TRNCYPXXISRXX',
    'TRNCYPXXLBNXX',
    'TRNCYPXXSYRXX',
    'TRNCYPXXTURXX',
    'TRNDEUXXGBRXX',
    'TRNDEUXXNORXX',
    'TRNDEUXXSWEXX',
    'TRNDJIXXYEMXX',
    'TRNDNKXXGBRXX',
    'TRNDNKXXNLDXX',
    'TRNDNKXXNORXX',
    'TRNDNKXXSWEXX',
    'TRNDZAXXESPXX',
    'TRNDZAXXFRAXX',
    'TRNDZAXXITAXX',
    'TRNEGYXXGRCXX',
    'TRNEGYXXJORXX',
    'TRNEGYXXSAUXX',
    'TRNERIXXSAUXX',
    'TRNERIXXYEMXX',
    'TRNESPXXMARXX',
    'TRNESTXXFINXX',
    'TRNFRAXXGBRXX',
    'TRNFRAXXIRLXX',
    'TRNGBRXXISLXX',
    'TRNGBRXXNLDXX',
    'TRNGBRXXNORXX',
    'TRNGRCXXITAXX',
    'TRNGRCXXLBYXX',
    'TRNGRLXXISLXX',
    'TRNIDNXXMYSXX',
    'TRNIDNXXSGPXX',
    'TRNINDSOLKAXX',
    'TRNINDWEOMNXX',
    'TRNINDWESAUXX',
    'TRNIRNXXOMNXX',
    'TRNITAXXMNEXX',
    'TRNITAXXTUNXX',
    'TRNJPNCEKORXX',
    'TRNJPNHORUSFE',
    'TRNJPNKYKORXX',
    'TRNLBYXXMLTXX',
    'TRNLTUXXSWEXX',
    'TRNLVAXXSWEXX',
    'TRNMARXXPRTXX',
    'TRNMLTXXTUNXX',
    'TRNMYSXXPHLXX',
    'TRNNLDXXNORXX',
    'TRNOMNXXPAKXX',
    'TRNPOLXXSWEXX',
    'TRNSAUXXSDNXX',
    'TRNSOMXXYEMXX',
    ]