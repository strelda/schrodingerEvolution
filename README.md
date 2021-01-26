# Časový vývoj jednorozměrného vlnového balíku

## Programy
**schrod.jl** (Julia, verze 1.5.3)

Módy:
- pohyb volné částice 
- testingMode - rozplývání Gaussovského balíku
- pohyb vlnového balíku v potenciálové jámě

Výsledek pro první dva případy je srovnáván s analytickým řešením. Na základě toho lze odhadnout parametry pro dostatečnou přesnost při výpočet třetího problému. Přesnost iterací je určena parametry
- časový krok *dt* 
- diskretizace prostoru na *n*-bodů. 
- škálování prostoru *scale*, aby *ψ[0]=ψ[n]=0*. 

Čas do kterého budou iterace prováděny je řízen fyzikálním časem *tExact*.


#### pomocné skripty
**imagesToVideo** - linux script pro převod obrázkové sekvence na video

**calculateItPlease** - linux script pro spuštění programu schrod.jl až po zobrazení výsledků


## Volná částice
Na následujícím obrázku vidíme počáteční stav vlnové funkce s parametry *scale=30*, *k=0.1*. ![img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_initial.jpeg?raw=true "vlastní stav")

Pro energii částice *En=5* a přesnost výpočtu danou *n=500*, *dt=5e-2* dostaneme v čase tExact=10 ![img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_500_30_5e-2_10_k0.1_E5.jpeg?raw=true "vlastní stav 1")

Časový vývoj je srovnán ve videích **videos/eigen_t=... .mp4** pro časové kroky *dt=t=1e-1* a *5e-2*. Pro *dt=1e-2* nebyl již rozdíl v daném čase *tExact* pozorovatelný.

Zvyšování přesnosti výpočtu lze ilustrovat na několika hodnotách, které program automaticky ukládá do souboru intDiff. Tato chyba je počítaná, jako <ψ-ψexact|ψ-ψexact>, ze kterého je navíc vynechaná rozumně velká oblast kolem krajů tak, aby braket nebyl ovlivněn okrajovými nepřesnostmi (viz spiky na předchozím obrázku, graf dole).
Tyto chyby se při konstantním *n* snižují, o 4 řády s každým řádem zmenšení *dt*. Tedy průměrná jednodimenzionální chyba *ψ-ψexact|<sub>x;</sub>* je skutečně úměrná *dt²*. 

#### Autokorelační funkce
*<ψ|ψ0>* je obecně komplexní funkce, která je v případě volné částice periodická. Perioda je v případě dt=1e-2 přibližně 23 časových kroků, tedy 23*dt=0.23 s. (kde uvažujeme sekundy, jako fyzikální jednotku času tExact)

V daném případě byla autokorelační funkce pro čas *tExact=10*
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_autocorrelation.jpeg?raw=true "Autokorelační funkce")


## Gaussovský balík
Vývoj Gaussovského balíku je silně ovlivněn okrajovými podmínkami a nepomůže ani jejich fixování na nulu (odraz vlny od pevného konce není nulový). Je tedy nutné nastavit dostatečně velkou škálu *scale*. 

Počáteční stav gaussovského balíku pro parametry *μ=1, ω=1.5, σ=0.5, p0=0, x0=-5* je
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_initial.jpeg?raw=true "Autokorelační funkce")


a po 3000 iteracích vyvine do stavu 
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_10000_20_1e-1_30_frame30.jpeg?raw=true "Autokorelační funkce")
kde lze opět pozorovat mírné odchýlení od přesného řešení, jelikož výpočet byl proveden s nízkými parametry *n=1000*, *dt=1e-1*.

V případě rozplývajícího se Gaussovského balíku není autokorelační funkce periodickou funkcí, ale vypadá
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelation.jpeg?raw=true "Autokorelační funkce")

#### Gaussovský balík v potenciálové jámě
Přejdeme-li však z testovacího módu (tedy *testingMode=1* v schrod.jl),čímž zapneme kvadratický potenciál, dostaneme např pro parametry *n=5000*, *dt=1e-2*, *scale=10* a čas *tExact=10* časový vývoj, viz **videos/gaussInPotential.mp4**. Autokorelační funkce pak bude
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelationPotential.jpeg?raw=true "Autokorelační funkce")

Perioda autokorelační funkce je v tomto případě přibližně 833 snímků (spočítáno na 24 periodách), tedy *Δt≈8.33* s.

###### Závislost autokorelační funkce na parametrech
Pro jiné hodnoty, než *μ=1*, např *μ=4*, můžeme pozorovat několik period modulovaných do jedné vlny
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelation_μ4.jpeg?raw=true "Autokorelační funkce")
avšak největší perioda funkce se nemění.

Parametru *ω* je frekvence autokorelační funkce úměrná, na parametru *σ* nezávisí, hybnost *p0* způsobí přebíhání částice v potenciálovém údolí, tedy kvůli přidané komplexitě pohybu periodu zvětší, viz **videos/pingpong.mp4** pro *p=1*. 