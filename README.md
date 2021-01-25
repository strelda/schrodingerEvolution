# Časový vývoj jednorozměrného vlnového balíku

hlavní program schrod.jl (Julia, verze 1.5.3)
imagesToVideo - linux script pro převod obrázkové sekvence ve video
calculateItPlease - linux script pro spuštění programu až po zobrazení výsledků

schrod.jl na požádání spočte časový vývoj volné částice, či rozplývání Gaussovského balíku. Výsledek je srovnáván s analytickým řešením, aby bylo možné určit chybu výpočtu. Ta závisí na časodém kroku *dt* a na diskretizaci prostoru na *n*-bodů. Dále je třeba škálovat prostor tak, aby se rozumná část grafu vešla na vyvíjenou oblast, k čemuž slouží parametr *scale*. Čas do kterého budou iterace prováděny je řízen fyzikálním časem *tExact*.


## Volná částice
Na obrázku images/eigen_initial.jpeg vidíme počáteční stav vlnové funkce s parametry *scale=30*, *k=0.1*. Pro energii částice En=5 a přesnost výpočtu danou n=500, dt=5e-5 která dostaneme v čase tExact=10 výsledek zobrazený na obrázku ![img](https://github.com/strelda/schrodingerEvolution/blob/main/images/eigen_500_30_5e-2_10_k0.1_E5.jpeg?raw=true "eigen").

Časový vývoj je srovnán ve videích schrodinger_t=... .mp4 pro časové kroky *dt=t=1e-1* a *5e-2*. Pro *dt=1e-2* nebyl již rozdíl v daném čase t*Exact* pozorovatelný.

Zvyšování přesnosti výpočtu lze ilustrovat na několika hodnotách, které program automaticky ukládá do souboru intDiff. Tato chyba je počítaná, jako <ψ-ψexact|ψ-ψexact>, ze kterého je navíc vynechaná rozumně velká oblast kolem krajů tak, aby braket nebyl ovlivněn okrajovými nepřesnostmi.
Tyto chyby se při konstantním *n* snižují, o 4 řády s každým řádem zmenšení *dt*. Tedy jednodimenzionální chyba *ψ-ψexact|<sub>x;</sub>* je skutečně úměrná *dt^2*. 

*<ψ|ψ0>* je obecně komplexní funkce, která je v případě volné částice periodická. Perioda je v případě dt=1e-2 přibližně 23 časových kroků, tedy 23*dt=1e-2
![Alt text](images/eigen_autoCorrelation.jpeg?raw=true "Autokorelační funkce")


##Gaussovský balík
Vývoj Gaussovského balíku je silně ovlivněn okrajovými podmínkami a nepomůže ani jejich fixování na nulu. Je tedy nutné nastavit dostatečně velkou škálu *scale*. 

Na obrázku images/gauss_initial.jpeg vidíme počáteční stav gaussovského balíku. Po čase *tExact=50* je tento stav vyvinut do stavu images/gauss_1000_60_1e-2_50.jpeg.

V případě Gaussovského balíku není autokorelační funkce periodickou funkcí, ale opět ji můžeme vykreslit, viz images/gauss_autoCorrelation.jpeg.

Zajímavé je pozorovat překmitávání vln na kraji Gaussovského balíku patrné např. z <ψ-ψexact|ψ-ψexact>.
