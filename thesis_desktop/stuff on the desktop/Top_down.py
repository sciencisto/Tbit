"""
Created on Fri Jan  3 14:29:00 2020
@author: kristine
"""
#Import modules
import math

def top_down_udskæring(udskæring, hals, højde, bredte):
    """Funtion der angiver opskriften på halsudskæringen
    
    - Udskæring:     valgmulighederne er her "thight" og "alm" hvor thight betyder at kraven vil være 
                     ved min hals mens alm betyder at kraven vil ligge omkring kravebenet, når alm og raglan
                     vælges anvendes der vendepinde. 
    - Hals:          valgmulighederne er "high" og "alm" hvor høj betyder at der strikkes længere i rib for
                     at der kan foldes uden syning, mens alm betyder at der strikkes en mindre ribkant som 
                     foldes."""
    #tætsiddende hals
    if udskæring == "thight":
        opslag = math.floor(bredte*4) # skal blive til 40 cm 
       
        if hals == "high":
            hals_længde = "14 cm"
            print("Tæt halsudskæring med rullekrave")
            print("---------------------------------")
            print("Slå %s masker op til udskæring og strik %s rib. Kraven skal blot foldes ud efter afsluttet stik, så juster denne længde efter hvordan du synes den skal ligge." % (opslag, hals_længde)) 
        
        else:
            hals_længde = "8 cm"
            print("Tæt Halsudskæring")
            print("-----------------")
            print("Slå %s masker op til udskæring og strik %s rib. Kraven skal foldes indad og sys, så den ligger dobbelt, hvorfor du med fordel kan justere denne længde så det passer dig." % (opslag, hals_længde)) 
    
    #udskæring med et mere løst fit
    elif udskæring == "alm":
        print("Almindelig Halsudskæring")
        print("-------------------------")
        opslag = math.floor(bredte*5.3) #skal blive til 53 cm 
        hals_længde = "6 cm"   
        print("Slå %s masker op til udskæring og strik %s rib. Kraven skal foldes indad og sys, så den ligger dobbelt, hvorfor du med fordel kan justere denne længde så det passer dig." % (opslag, hals_længde))    
        
    else:
        print("Den halsudskæring kan endnu ikke vælges")

def top_down_bærestykke(udskæring, bærestykke, højde, bredte):
    """Funktion der udskiver opskriften på bærestykket.
    Valgmulighederne er "raglan" og "rundt" hvor raglan betyder at der laves karakteristiske
    udtagninger ved skuldre. rundt betyder, at der laves løbende udtagninger og der dermed ikke
    opstår de karakteristiske udtagninger."""
    
    if udskæring == "thight":
        opslag = math.floor(bredte*4)
        if bærestykke == "raglan":
            skulder = math.floor(opslag/6)
            front = math.ceil(opslag/6 * 2)
            masker_ved_afslut = bredte*16
            
            print("Raglan Bærestykke")
            print("-----------------")
            print("Sæt en markør ved omganges start, dernæst efter %s masker. Dette angiver hvor den ene skulder vil være. Dernæst fortsættes glaststrik %s masker hvilket nu er fronten af sweateren. %s masker afsættes til den anden skulder og der burde være %s masker tilbage til bagsiden før omgangen er slut. Der er dermed stadig %s masker på pinden. Nu startes bærestykket. Der udtages raglanmasker ved de fire markøre på hveranden rundte indtil raglan maske linjen er minumum 22 cm fra glatstrikens start. Der bør være ca. %s masker ved aflutningen af bærestykket" % (skulder, front, skulder, front, opslag, masker_ved_afslut))
    
        elif bærestykke == "dyb raglan":
            skulder = math.floor(opslag/6)
            front = math.ceil(opslag/6 * 2)
            masker_ved_afslut = bredte*16+(2*(math.ceil((bredte/10)*4)))
            
            print("Dyb Raglan Bærestykke")
            print("-----------------")
            print("Sæt en markør ved omganges start, dernæst efter %s masker. Dette angiver hvor den ene skulder vil være. Dernæst fortsættes glaststrik %s masker hvilket nu er fronten af sweateren. %s masker afsættes til den anden skulder og der burde være %s masker tilbage til bagsiden før omgangen er slut. Der er dermed stadig %s masker på pinden. Nu startes bærestykket. Der udtages raglanmasker ved de fire markøre på hveranden rundte indtil der er ca. %s masker ved aflutningen af bærestykket, omkring 25 cm" % (skulder, front, skulder, front, opslag, masker_ved_afslut))
            
        elif bærestykke == "rundt": 
            print("Rundt Bærestykke")
            print("-----------------")
            udtagsantal = math.ceil(((højde/10)*3.5))*4
            masker_ved_afslut = bredte*16
            
            print("Sæt en markør ved omganges start, dernæst laves den første pind i glatstrik. På anden pind laves samlet %s udtagninger. Dette gøres ved at tilføje en maske efter hver anden maske, sæt endnu en markør under omgangsmarkøreren for at angive hvornår udtagningerne skete. Strik 3.5 cm. Nu tilføjes en maske efter hver 3. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 4. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 5. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 6. maske. Sæt en markør. strik 3.5 cm og tilføj en maske efter hver 7. maske. Strik til bærestykket er 22 cm. juster til %s masker om nødvendigt" % (udtagsantal, masker_ved_afslut))

        elif bærestykke == "sunday": 
            print("Sunday Bærestykke")
            print("-----------------")
            udtagsantal = math.ceil(((højde/10)*3.5))*4
            masker_ved_afslut = bredte*16
            
            print("Sæt en markør ved omganges start, dernæst laves den første pind i glatstrik. På anden pind laves udtagninger således at 1 ret en vrang mønsteret nu bliver 2 ret 1 vrang. Dette gøres ved at tilføje en maske efter hver anden maske, sæt endnu en markør under omgangsmarkøreren for at angive hvornår udtagningerne skete. Strik 3 cm. Nu tilføjes en maske efter hver vrang-maske så mønsteret bliver 2 ret 2 vrang. Sæt en markør. Strik 4 cm. Nu tilføjes en maske mellem de to retmasker. Sæt en markør, mønsteret er nu 3 ret og 2 vrang. Strik 4 cm. Nu tilføjes en maske mellem hver vrangmaskerne, sådan mønsteret bliver 3 ret, 3 vrang. Sæt en markør. Strik 5 cm. Nu tilføjes enndu en ret-maske, sæt en markør. Mønsteret er nu 4 ret og 3 vrang. strik 5 cm. Tilføj den sidste ekstra maske, endnu en vrang så mønsteret er 4 ret, 4 vrang, strik til bærestykket er 26 cm fra kraven")

    elif udskæring == "alm":
        opslag = math.floor(bredte*5.3)
        if bærestykke == "raglan":
            skulder = math.floor(opslag/6)
            front = math.ceil(opslag/6 * 2)
            masker_ved_afslut = bredte*16
            
            print("Raglan Bærestykke")
            print("-----------------")
            print("Sæt en markør ved omganges start, dernæst efter %s masker. Dette angiver hvor den ene skulder vil være. Dernæst fortsættes glaststrik %s masker hvilket nu er fronten af sweateren. %s masker afsættes til den anden skulder og der burde være %s masker tilbage til bagsiden før omgangen er slut. Der er dermed stadig %s masker på pinden. Nu startes bærestykket. Der udtages raglanmasker ved de fire markøre på hveranden rundte i 2 cm. Dernæst startes vendepinde, der laves %s vendepinde (german short rows). Vendepinde og raglanudtagninger: første vendepind går 2 masker ind på forstykket, den næste 4 og så vidre. Raglanudtagninger foregår på hver anden som sædvantigt. Derefter strikkes der rundt indtil raglan maske linjen er minumum 22 cm fra glatstrikens start. Der bør være ca. %s masker ved aflutningen af bærestykket" % (skulder, front, skulder, front, opslag, math.ceil(højde/2), masker_ved_afslut))
    
        elif bærestykke == "dyb raglan":
            skulder = math.floor(opslag/6)
            front = math.ceil(opslag/6 * 2)
            masker_ved_afslut = bredte*16+(2*(math.ceil((bredte/10)*4)))
            
            print("Dyb Raglan Bærestykke")
            print("-----------------")
            print("Sæt en markør ved omganges start, dernæst efter %s masker. Dette angiver hvor den ene skulder vil være. Dernæst fortsættes glaststrik %s masker hvilket nu er fronten af sweateren. %s masker afsættes til den anden skulder og der burde være %s masker tilbage til bagsiden før omgangen er slut. Der er dermed stadig %s masker på pinden. Nu startes bærestykket. Der udtages raglanmasker ved de fire markøre på hveranden rundte i 2 cm. Dernæst startes vendepinde, der laves %s vendepinde (german short rows). Vendepinde og raglanudtagninger: første vendepind går 2 masker ind på forstykket, den næste 4 og så vidre. Raglanudtagninger foregår på hver anden som sædvantigt. Derefter strikkes der rundt indtil raglan maske linjen er minumum 25 cm fra glatstrikens start. Der bør være ca. %s masker ved aflutningen af bærestykket" % (skulder, front, skulder, front, opslag, math.ceil(højde/2), masker_ved_afslut))
            
        elif bærestykke == "rundt": 
            print("Rundt Bærestykke")
            print("-----------------")
            udtagsantal = math.ceil(((højde/10)*3.5))*4
            masker_ved_afslut = bredte*16
            print("Sæt en markør ved omganges start, dernæst laves den første pind i glatstrik. På anden pind laves samlet %s udtagninger. Dette gøres ved at tilføje en maske efter hver anden maske, sæt endnu en markør under omgangsmarkøreren for at angive hvornår udtagningerne skete. Nu laves der vendepinde, german shortrows på ryg siden, dette gøres på de næste 3.5 cm indtil næste udtagning. første vendepind går 2 masker ind på forstykket, den næste 4 og så vidre. Nu tilføjes en maske efter hver 3. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 4. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 5. maske. Sæt en markør. Strik 3.5 cm. Nu tilføjes en maske efter hver 6. maske. Sæt en markør. strik 3.5 cm og tilføj en maske efter hver 7. maske. Strik til bærestykket er 22 cm. juster til %s masker om nødvendigt" % (udtagsantal, masker_ved_afslut))

        elif bærestykke == "sunday": 
            print("Sunday Bærestykke")
            print("-----------------")
            udtagsantal = math.ceil(((højde/10)*3.5))*4
            masker_ved_afslut = bredte*16
            
            print("Sæt en markør ved omganges start, dernæst laves den første pind i glatstrik. På anden pind laves udtagninger således at 1 ret en vrang mønsteret nu bliver 2 ret 1 vrang. Dette gøres ved at tilføje en maske efter hver anden maske, sæt endnu en markør under omgangsmarkøreren for at angive hvornår udtagningerne skete. Strik 3 cm. Nu tilføjes en maske efter hver vrang-maske så mønsteret bliver 2 ret 2 vrang. Sæt en markør. Strik 4 cm. Nu tilføjes en maske mellem de to retmasker. Sæt en markør, mønsteret er nu 3 ret og 2 vrang. Strik 4 cm. Nu tilføjes en maske mellem hver vrangmaskerne, sådan mønsteret bliver 3 ret, 3 vrang. Sæt en markør. Strik 5 cm. Nu tilføjes enndu en ret-maske, sæt en markør. Mønsteret er nu 4 ret og 3 vrang. strik 5 cm. Tilføj den sidste ekstra maske, endnu en vrang så mønsteret er 4 ret, 4 vrang, strik til bærestykket er 26 cm fra kraven på rygsiden")


    else: 
        print("der gik noget galt i print af bærestykket")
        
def top_down_krop(bærestykke,højde, bredte):
    armhule_opslag = math.ceil((bredte/10)*5)
    masker_ved_afslut = bredte*16
        
    print("Krop")
    print("------")
    print("Nu startes kroppen. Dette gøres ved at sætte ærmer til hvile og lave nye masker under armen.")
    
    if bærestykke == "raglan":
        print("Sæt maskerne i skulder sektionerne til hvile. Inkluder ikke selve raglan masken. Slå %s masker op under armene" % (armhule_opslag))
        print("Der stikkes i glatstrik indtil kroppen måler 35 cm fra armhulen, og der afsluttes med at strikkes en ribkant, 4 cm. Vurder om der skal laves indtagninger i denne. Luk forsigtigt af.")

    
    elif bærestykke == "dyb raglan":
        print("Sæt maskerne i skulder sektionerne til hvile. Inkluder ikke selve raglan masken. undlad at slå masker op under armen")
        print("Der stikkes i glatstrik indtil kroppen måler 35 cm fra armhulen, og der afsluttes med at strikkes en ribkant, 4 cm. Vurder om der skal laves indtagninger i denne. Luk forsigtigt af.")
    
    elif bærestykke == "rundt":
        ærmer = math.floor(masker_ved_afslut/6)
        krop = math.ceil((masker_ved_afslut/6)*2)
        print("Sæt %s masker til hvile som det ene ærme og strik %s masker frem til der kan sættes %s masker til hvile. Her slås %s masker op i armhulen før ryggen strikkes." % (ærmer, krop, ærmer, armhule_opslag))
        print("Der stikkes i glatstrik indtil kroppen måler 35 cm fra armhulen, og der afsluttes med at strikkes en ribkant, 4 cm. Vurder om der skal laves indtagninger i denne. Luk forsigtigt af.")

    elif bærestykke == "sunday":
        ærmer = math.floor(masker_ved_afslut/6)
        krop = math.ceil((masker_ved_afslut/6)*2)
        print("Sæt %s masker til hvile som det ene ærme og strik %s masker frem til der kan sættes %s masker til hvile. Her slås %s masker op i armhulen før ryggen strikkes." % (ærmer, krop, ærmer, armhule_opslag))
        print("Der stikkes i enten glatstrik eller 4 ret 4 vrang mønsteret fortsætter indtil kroppen måler 35 cm fra armhulen, og der afsluttes med at strikkes en ribkant, 4 cm. Vurder om der skal laves indtagninger i denne. Luk forsigtigt af.")

    else: 
        print("Der er sket en fejl i indtastningen, prøv igen")
        
    
def top_down_ærmer(ærmer, højde, bredte):
    armhule_opsamling = 2 + (math.ceil((bredte/10)*4))
    print("Ærmer")
    print("------")
    
    if ærmer == "thight":
        print("Der samles %s masker op i armhulen og der strikkes i glat eller mønstret strik efter behov rundt*. Husk at strik to masker sammen i hver side af åbningen for at undgå huller. Hver femte centimeter strikkes to masker samen i begyndelsen og slutningen af rundten under armen. Fortsæt indtil ærmet måler minumum 44 cm, strik en ribkant med %s masker, 4 cm. Luk forsigtigt af." % (armhule_opsamling, bredte*2))
    
    elif ærmer == "loose":
        print("Der samles %s masker op i armhulen og der strikkes i glat eller mønstret strik efter behov rundt*. Husk at udvide strikke to masker sammen i hver side af åbningen for at undgå huller. Strik 10 cm lav den første udtagningmed to masker i begyndelsen og slutningen af rundten under armen, ligesom ved raglan udtagninger. Dernæst strikkes 5 cm. Endnu en udtagning. yderligere 5 cm enndu en udtagning. Dette er den sidste. Strik 20 cm glatstrik. Ærmet måler nu 40 cm. Der laves indtagninger på 3 runder fordelt hver %s pind, hvor 2 masker strikkes sammen i starten og afslutningen af omgangen på samme linje som den orginale udtagning. Stik parvis maskerne sammen og afslut med en ribkant med %s masker, 4 cm. Luk forsigtigt af." % (armhule_opsamling, math.ceil((højde/2)/3), bredte*2))        
        
    else:
        print("Den type ærmer kan endnu ikke vælges")
        
    print("* Noter at mønstret strik f.eks. kan være perlestrik, rib, hulmønster eller lignende.")

def garn_beregner(løbelængde, pind):
    if pind in range(2,4):
        nøgler = math.ceil(1000/løbelængde)
    elif pind in range(4,6):
        nøgler = math.ceil(800/løbelængde)
    elif pind in range(6,8):
        nøgler = math.ceil(650/løbelængde)
    elif pind in range(8,10):
        nøgler = math.ceil(400/løbelængde)
    elif pind in range(10,12):
        nøgler = math.ceil(300/løbelængde)
    return nøgler

    
def top_down_sweater(udskæring, hals, bærestykke, ærmer, højde, bredte, pind, løbelængde):
    print("Du har bedt om opskiften på en sweater i Kristine størrelse med følgende specifikationer: %s udskæring med en %s hals, et %s bærestykke samt %s ærmer. Din strikkeprøve målte %s pinde og %s masker på 10*10 cm" % (udskæring, hals, bærestykke, ærmer, højde, bredte))
    
    print("---------------------------------------------")
    
    print("Opskriften er som følger")   
    top_down_udskæring(udskæring, hals, højde, bredte)   
    
    top_down_bærestykke(udskæring, bærestykke, højde, bredte)   
    
    top_down_krop(bærestykke,højde, bredte)   
    
    top_down_ærmer(ærmer, højde, bredte)
   
    print("Dette er opskiften, det er esitmeret at du skal anvende den følgende mængde garn (antal nøgler):")
    nøgler = garn_beregner(løbelængde, pind)
    return nøgler

"""

    
top_down_sweater(udskæring="thight", hals="alm", bærestykke="raglan", ærmer="loose", højde = 14, bredte = 11)
Tæt Halsudskæring
-----------------
Slå 44 masker op til udskæring og strik 8 cm rib. Skift til glatstrik og sæt følgende markøre på den første rundte:
Raglan Bærestykke
-----------------
ved omganges start, dernæst efter 7 masker. Dette angiver hvor den ene skulder vil være. Dernæst fortsættes glaststrik 15 masker hvilket nu er fronten af sweateren. 7 masker afsættes til den anden skulder og der burde være 15 masker tilbage til bagsiden før omgangen er slut. Der er dermed stadig 44 masker på pinden. Nu startes bærestykket. Der udtages raglanmasker ved de fire markøre på hveranden rundte indtil raglan maske linjen er 22 cm fra glatstrikens start. Der bør være ca. 176 masker ved aflutningen af bærestykket
Krop
------
Nu startes kroppen. Dette gøres ved at sætte ærmer til hvile og lave nye masker under armen.
Sæt maskerne i skulder sektionerne til hvile. inkluder ikke selve raglan masken. Slå 5 masker op under armene
Der stikkes i glatstrik indtil kroppen måler 28 cm fra armhulen, og der afsluttes med at strikkes en ribkant, 4 cm. Vurder om der skal laves indtagninger i denne. Luk forsigtigt af.
Ærmer
------
Der samles et passende antal masker op i armhulen og der strikkes i glat eller mønstret strik efter behov rundt. Husk at strik to masker sammen i hver side af åbningen for at undgå huller. Hver femte centimeter øges maskeabtalet med to masker samen i begyndelsen og slutningen af rundten under armen. Fortsæt indtil ærmet måler minumum 40 cm, lav dernæst indtagninger og afslut med en ribkant med 22 masker, 4 cm. Luk forsigtigt af.
Dette er opskriften på din sweater /Kristine       
""" 

        
