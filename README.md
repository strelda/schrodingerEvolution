# Časový vývoj jednorozměrného vlnového balíku

hlavní program schrod.jl (Julia, verze 1.5.3)
imagesToVideo - linux script pro převod obrázkové sekvence ve video
calculateItPlease - linux script pro spuštění programu až po zobrazení výsledků

schrod.jl na požádání spočte časový vývoj volné částice, či rozplývání Gaussovského balíku. Výsledek je srovnáván s analytickým řešením, aby bylo možné určit chybu výpočtu. Ta závisí na časodém kroku *dt* a na diskretizaci prostoru na *n*-bodů. Dále je třeba škálovat prostor tak, aby se rozumná část grafu vešla na vyvíjenou oblast, k čemuž slouží parametr *scale*. Čas do kterého budou iterace prováděny je řízen fyzikálním časem *tExact*.


## Volná částice
Na následujícím obrázku vidíme počáteční stav vlnové funkce s parametry *scale=30*, *k=0.1*. ![img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_initial.jpeg?raw=true "vlastní stav")

Pro energii částice *En=5* a přesnost výpočtu danou *n=500*, *dt=5e-5* dostaneme v čase tExact=10 ![img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_500_30_5e-2_10_k0.1_E5.jpeg?raw=true "vlastní stav 1").

Časový vývoj je srovnán ve videích **videos/eigen_t=... .mp4** pro časové kroky *dt=t=1e-1* a *5e-2*. Pro *dt=1e-2* nebyl již rozdíl v daném čase *tExact* pozorovatelný.

Zvyšování přesnosti výpočtu lze ilustrovat na několika hodnotách, které program automaticky ukládá do souboru intDiff. Tato chyba je počítaná, jako <ψ-ψexact|ψ-ψexact>, ze kterého je navíc vynechaná rozumně velká oblast kolem krajů tak, aby braket nebyl ovlivněn okrajovými nepřesnostmi.
Tyto chyby se při konstantním *n* snižují, o 4 řády s každým řádem zmenšení *dt*. Tedy jednodimenzionální chyba *ψ-ψexact|<sub>x;</sub>* je skutečně úměrná *dt^2*. 

*<ψ|ψ0>* je obecně komplexní funkce, která je v případě volné částice periodická. Perioda je v případě dt=1e-2 přibližně 23 časových kroků, tedy 23*dt=0.23. 

Autokorelační funkce:
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_autocorrelation.jpeg?raw=true "Autokorelační funkce")


##Gaussovský balík
Vývoj Gaussovského balíku je silně ovlivněn okrajovými podmínkami a nepomůže ani jejich fixování na nulu. Je tedy nutné nastavit dostatečně velkou škálu *scale*. 

Počáteční stav gaussovského balíku
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_initial.jpeg?raw=true "Autokorelační funkce")


se po čase 30 iteracích vyvine do stavu 
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_10000_20_1e-1_30_frame30.jpeg?raw=true "Autokorelační funkce"),
kde lze opět pozorovat mírné odchýlení od přesného řešení, jelikož výpočet byl proveden s nízkými parametry *n=1000*, *dt=1e-1*.

V případě rozplývajícího se Gaussovského balíku není autokorelační funkce periodickou funkcí, ale vypadá
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelation.jpeg?raw=true "Autokorelační funkce")

Přejdeme-li však z testovacího módu (tedy *testingMode=1* v schrod.jl),čímž zapneme kvadratický potenciál, dostaneme např pro parametry *n=5000*, *dt=1e-2*, *scale=10* a čas *tExact=10* časový vývoj, viz **videos/gaussInPotential.mp4**. Autokorelační funkce pak bude
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelationPotential.jpeg?raw=true "Autokorelační funkce")

Perioda autokorelační funkce je v tomto případě přibližně 833 snímků, tedy *Δt≈8.33* s.

####Závislost autokorelační funkce na parametrech
Pro jiné hodnoty, než *μ=1*, např *μ=4* viz obrázek můžeme pozorovat několik period modulovaných do jedné vlny
![Img](https://github.com/strelda/schrodingerEvolution/blob/main/images/gauss_autocorrelation_μ4.jpeg?raw=true "Autokorelační funkce")
avšak největší perioda funkce se nemění.

Parametru *ω* je frekvence úměrná, na parametru *σ* nezávisí, hybnost *p* způsobí přebíhání částice v potenciálovém údolí, tedy kvůli přidané komplexitě pohybu periodu zvětší, viz **videos/pingpong.mp4** pro *p=1*. 