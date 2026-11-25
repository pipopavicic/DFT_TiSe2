# DFT fononi monosloja 1T-TiSe₂

Otvori **Vodic_DFT_TiSe2.pdf** za vezu s vježbama 2, 4, 5 i 6. Dvije stranice s naslovima „Što si zbrojio u bilježnici?” i „Kako usporediti dvije amplitude?” polaze od tvoje tablice pomaka Ti i spremljenih MACE minimuma. Zadnja stranica vodi kroz stvarne naredbe. DFT dio koristi Quantum ESPRESSO i PBE za sve temperature.

## Od jednog M-moda do njihova zbroja

U dijelu 9 bilježnice crtaš energiju prema amplitudi pomaka. Jednim smjerom pomičeš atome prema M₁, a drugim prema zbroju M₁ + M₂ + M₃. Oba puta kreću od iste visokosimetrične strukture. Fononska disperzija prikazuje frekvencije prema valnom vektoru i za tu je početnu strukturu zajednička. U bilježnici oba izvoza nose oznaku `M1`, pa drugi može prepisati prvi. U DFT mapi svaki put dobiva vlastitu datoteku.

| Spremljeni rezultat iz bilježnice | M₁, dubina (meV/2×2) | 3Q, dubina (meV/2×2) | Iznos najvećeg pomaka atoma u minimumu M₁ / 3Q (Å) |
|---|---:|---:|---:|
| MACE-PBE | 11,1886 | 19,3714 | 0,056 / 0,084 |
| MACE-r²SCAN | 9,6650 | 16,5313 | 0,048 / 0,072 |

Ovo su očitanja minimuma na mreži MACE pomaka. DFT krivulje za 100, 200, 400, 800, 1200 i 1600 K nastaju kad pokreneš račune. U svakoj točki dobivaš razliku elektronskih slobodnih energija ΔF = F(A,T) − F(0,T).

## Rezultati na ovom računalu

Za osnovne postavke otvori **racuni/k24_q6_17503ad3bf/rezultati**. Datoteka `status.csv` pokazuje do koje se faze stiglo. `usporedba.csv` dobiva broj tek kada odgovarajući račun završi. Za prvi dovršeni ionski račun tu je i `relaksacija_1x1.csv`. Završni položaji, izlazi QE-a i svojstveni vektori ostaju u mapi pojedine temperature. Mapa s oznakom **k4_q2** pripada tehničkoj provjeri pri 1600 K na gruboj mreži.

Primjer čitanja stvarnog izlaza nalazi se u `racuni/k4_q2_b4968b2f3c/T1600K/kontrola_Gamma_1x1/frekvencije.csv`. U Γ primitivne ćelije nakon primjene akustičkog pravila sume dobiveno je 0, 0, 0, 3,924, 3,924, 4,971, 4,971, 5,696 i 8,582 THz. Tri nule odgovaraju zajedničkim translacijama svih atoma. U ovoj ćeliji treba otići do M da bi se ispitao uzorak koji mijenja predznak između susjednih ćelija.

Na istoj probnoj disperziji minimum je u M. Izravni zapis za simetrijski ekvivalentnu M-točku fononske mreže daje −1,087780 THz, a interpolacija duž putanje −1,087779 THz. Slaganje provjerava putanju, jedinice i interpolaciju. Za fizikalno tumačenje pri 1600 K valja usporediti glavnu i gušće elektronske mreže. Graf je u `racuni/k4_q2_b4968b2f3c/rezultati/disperzija_1600K.png`.

## Ažuriranje koda Gitom na Srcu

Repozitorij sadrži skripte, postavke, pseudopotencijale i vodič. Radne datoteke i izračunati rezultati u `racuni/` ostaju na stroju na kojem pokrećeš Quantum ESPRESSO. Nove verzije koda preuzmi između pokretanja studije.

Ako već imaš mapu `~/zadatak/DFT_TiSe2` s rezultatima, možeš je jednom povezati s repozitorijem. Pripremi adresu repozitorija koji sadrži ovaj commit na grani `main`. Ako koristiš prijenosnu datoteku, učitaj `DFT_TiSe2.bundle` u `~/zadatak/`. U terminalu pokreni donji blok. Kada se pojavi pitanje o izvoru, za učitani bundle pritisni Enter ili zalijepi adresu svojeg repozitorija.

```bash
(
    set -eu
    cd ~/zadatak/DFT_TiSe2
    if git rev-parse --git-dir >/dev/null 2>&1; then
        echo "Ova je mapa već u Git repozitoriju. Provjeri git status i git remote -v."
        exit 1
    fi
    read -r -p "URL repozitorija ili Enter za DFT_TiSe2.bundle: " izvor
    izvor=${izvor:-../DFT_TiSe2.bundle}
    git ls-remote --exit-code "$izvor" refs/heads/main >/dev/null
    git init -b main
    git remote add origin "$izvor"
    git fetch origin main
    git reset --mixed FETCH_HEAD
    git -c user.name=Jupyter -c user.email=jupyter@localhost stash push -m "Lokalne skripte prije povezivanja"
    git branch --set-upstream-to=origin/main main
    git status --short --branch
)
```

Ako su se tvoje skripte ili `postavke.json` razlikovali od verzije u Gitu, Git sprema te razlike u `stash`. Pogledaj `git stash list`, a zatim `git stash show -p` ako je stvoren novi zapis. Prije sljedećeg računa provjeri vlastite izmjene postavki jer one određuju mapu rezultata. Osobne bilježnice i postojeće mape izračuna ostaju gdje su bile.

Kad `origin` pokazuje na udaljeni repozitorij, nove verzije preuzimaš naredbom `git pull --ff-only`. Vlastite izmjene šalješ nakon `git add ime_datoteke` i `git commit` naredbom `git push`. Ako si počeo iz bundlea, jednom postavi `git remote set-url origin URL_REPOZITORIJA`. Za novi prazan udaljeni repozitorij prva je objava `git push -u origin main`.

## Pokretanje

Iz mape projekta MRMFKT na ovom Windows računalu:

```powershell
cd DFT_TiSe2
..\.venv\Scripts\python.exe studija.py status
..\.venv\Scripts\python.exe studija.py run
```

`run` automatski pripremi ulaze, koristi dostupne potvrđene korake i obrađuje svih šest temperatura. Nakon svake obrađene temperature ili dijela serije osvježava tablice i grafove. Ako želiš ponovno generirati prikaze iz dosad dovršenih računa, pokreni `..\.venv\Scripts\python.exe studija.py report`. Dvije instance cijele serije za iste postavke zaštićene su zaključavanjem. Tijekom jednog dugog QE koraka napredak vidiš u njegovu izlazu, primjerice `T0100K/fononi_1x1/M.out`. Ako terminal ili računalo prekine račun, sljedeći `run` kreće od zadnjeg potvrđenog koraka. Prekinuti QE korak ponavlja se.

U standardnom korisničkom PowerShellu cijelu seriju možeš pokrenuti i s `.\Pokreni.ps1`. Opcija `-Pozadina` otvara proces bez prozora i zapisuje dnevnik `serija_*.log`. Na udaljenom računalu koristi svoj uobičajeni red poslova za dugotrajne račune. Program radi na CPU-u. Instalacija MACE-a s CUDA-om odnosi se na prvi dio projekta.

Za praćenje već pokrenute cijele serije iz drugog terminala postoji i `python studija.py watch`. Osvježava `PREGLED.md`, tablice i grafove čim pojedina faza završi i izlazi kada se serija zaustavi.

### Ako prenosiš paket na drugo računalo

Paket sadrži oba točna i provjerena UPF pseudopotencijala. U novom Python okruženju instaliraj ovisnosti s `python -m pip install -r requirements.txt`. Na Windowsu naredba `python instaliraj_qe_windows.py` preuzima provjereno prijenosno izdanje QE 7.5 u lokalnu mapu `engine`. Na Linuxu je dovoljan Quantum ESPRESSO s programima `pw.x`, `ph.x`, `q2r.x`, `matdyn.x` i `dynmat.x` u putanji. Varijabla `QE_BIN` može označiti i posebnu mapu tih programa. Na ovom računalu prijenosno izdanje već je pripremljeno.

## Kako proći gradivo postupno

Ove naredbe koriste istu konfiguraciju i iste podatke kao cijela serija. U primjerima je `python` interpreter projektnog okruženja.

| Što želiš proučiti | Naredba unutar ove mape |
|---|---|
| Vidjeti sve polazne QE ulaze i temperature | `python studija.py prepare` |
| Ionska relaksacija visokosimetrične ćelije pri 100 K | `python studija.py run --temperatures 100 --part relax1` |
| Devet frekvencija izravno u M | `python studija.py run --temperatures 100 --part m` |
| Cijela disperzija Γ–M–K–Γ | `python studija.py run --temperatures 100 --part mesh` |
| 36 svojstvenih parova u Γ simetrične 2×2 i vektori M₁, M₂, M₃, 3Q | `python studija.py run --temperatures 100 --part ref` |
| Dva presjeka, svaki na 13 amplituda s jednom zajedničkom nulom | `python studija.py run --temperatures 100 --part pes` |
| Slobodna 2×2 relaksacija i novi Γ fononi | `python studija.py run --temperatures 100 --part rel` |
| Provjeriti savijanje, normiranje zbroja, zakrivljenosti i QE formate na analitičkom primjeru | `python provjera.py` |
| Izvršiti testove konvergencije za M | `python konvergencija.py --run --temperatures 100 1600` |

Na najnižoj temperaturi izdvajamo najnižu granu iz svakog od triju M-sektora. Predznak biramo kao u bilježnici prema najvećoj apsolutnoj komponenti u prvoj primitivnoj kopiji. Na idućoj temperaturi u istom sektoru tražimo najveći apsolutni preklop i usklađujemo predznak. Tako možemo pratiti iste grane i kada postanu stabilne. Izdvajanje smjerova samo na višoj temperaturi može prvo izračunati nedostajuće prethodne reference. Ako je na početku M₁ stabilan i drugi sektor sadrži nestabilnost, kod dodatno radi presjek duž najnestabilnijeg sektora. Svi nestabilni vektori pune dijagonalizacije dostupni su u `nestabilni_vektori.csv`.

### Kako čitati dva presjeka

Za prvi pregled rezultata `--part ref` otvori **pomaci_M1_0100K.png** i **pomaci_3Q_0100K.png** u podmapi `rezultati`. Svaka slika prikazuje položaje svih 12 atoma sa strelicama i istim mjerilom uvećanja. Komponentu izvan ravnine možeš očitati u CSV-u. Tablica `M_triplet_prema_T.csv` prati frekvencije tri M-grane, raspon njihovih svojstvenih vrijednosti, 3Q zakrivljenost i rezidual kroz temperature.

- U `gamma_2x2_ref/smjerovi_M.json` za svaki M-mod nalaze se predznačena frekvencija u THz i svojstvena vrijednost u predznačenim THz². Kombinacija ima jednake koeficijente 1/√3 za tri maseno ponderirana, ortonormirana vektora. Zapisujemo projiciranu zakrivljenost eᵀDe i rezidual ‖De − (eᵀDe)e‖. Za točno degeneriran triplet i kombinacija je svojstveni vektor. Pri razdvojenim vrijednostima presjek i dalje ima dobro definiranu zakrivljenost.
- U `smjer_M1.csv`, `smjer_M2.csv`, `smjer_M3.csv` i `smjer_3Q.csv` stupci `e_*` daju maseno ponderirane komponente. Stupci `u_po_Q_*` dijele ih s √m atoma, a `u_*_rel` skaliraju taj pomak na najveću atomsku duljinu jedan. U DFT CSV-u `A_A` je potpisani najveći atomski pomak u Å. `Q_sqrt_amu_A` je jedinično normirana masena koordinata. Za isti fizički pomak u tvojoj ćeliji s nenormiranim zbrojem vrijedi Q_DFT-konvencija = √3 × Q_bilježnica. Na toj izvornoj osi harmonijska parabola kombinacije treba zbroj triju svojstvenih vrijednosti; na normiranoj osi treba njihov prosjek kada su translacijski sektori razdvojeni.
- Svaki presjek ima zaseban `pes/M1/potencijal.csv` ili `pes/3Q/potencijal.csv` sa stvarnom ΔF i harmonijskom procjenom. Oba dijele točno isti F(0,T). Zadani koraci idu od −0,16 do +0,16 Å, što za ova dva puta znači 24 nenulta SCF računa po temperaturi. Obje strane svakog presjeka računaju se neovisno. U `zakrivljenost.json` usporedi središnju razliku iz energija s DFPT projekcijom, u eV/Å² po 2×2. Uspoređuj istu apscisu i istu jedinicu energije.
- U `rezultati/PES_M1_3Q_0100K.png` nalaze se oba presjeka za 100 K. `PES_usporedba.csv` sadrži dubinu na uzorkovanoj mreži, amplitudu najniže točke, unutarnje minimume, rubnu oznaku i asimetriju. Nula kao najniža točka daje dubinu nula. Najniža vrijednost na rubu znači da položaj minimuma treba dodatno uzorkovati. Mreža DFT amplituda rjeđa je od one u bilježnici (61 MACE točka između ±0,12 Å).
- Slobodna 2×2 relaksacija na svakoj temperaturi kreće od male slučajne distorzije. `rezultati/distorzija_M.csv` daje tri predznačene projekcije konačne statičke distorzije (Å√amu) kada je dostupna i simetrična referenca. `grane_relaksirane.csv` prikazuje uz frekvenciju svake grane i projekcije svojstvenog vektora na odabrani 3Q smjer i na cijeli M-triplet. Kopije su između temperatura usklađene primitivnom translacijom; referentni M-vektori imaju fiksno ishodište. Projekcija na jedan 3Q smjer ovisi o domeni. `najveci_pomak.json` i `relaksirani_modovi.png` izdvajaju praćenu granu najvećeg raspona frekvencije nakon najmanje dvije obrađene temperature.

## Rad s ograničenim diskovnim prostorom

Na Jupyteru se u zasebnoj ćeliji može uključiti opcija za manju potrošnju diska. Za cijelu studiju pokreni

```bash
%%bash
set -e
cd ~/zadatak/DFT_TiSe2
export TISE2_LOW_DISK=1
python -u studija.py run
```

Uz zadanog jednog radnika skripta dovrši M, disperziju, potencijale i relaksiranu strukturu pojedine temperature prije prijelaza na sljedeću. Nakon pogreške staje kako bi se mogao pregledati izlaz. Ako već imaš ćeliju s petljama i odabranim fazama (`--part`), dodaj isti `export` nakon `cd` i zadrži svoj redoslijed. Ako je disk već pun, prvo oslobodi prostor. Promjena vrijedi za novo pokretanje programa.

`ph.x` dobiva `reduce_io=.true.`. Za relaksacije i pojedinačne energije potencijala `pw.x` dobiva `disk_io='nowf'`. SCF iz kojeg slijede fononi dobiva `disk_io='low'` i sprema valne funkcije koje će čitati `ph.x`. Te opcije mogu povećati potrošnju RAM-a. Prekinuti fononski korak s `reduce_io` treba pokrenuti ispočetka. Skripta prepoznaje dovršene korake iz oba načina pohrane. Stvarne postavke pojedinog pokretanja vidiš u datoteci `.in` uz njegov izlaz.

Periodični ispis u bilježnici i `rezultati/status.log` prikazuju i trenutačno slobodan prostor u GiB. Taj ispis uzorkuje stanje svakih 20 minuta. Za provjeru u bilo kojem trenutku u terminalu upiši `df -h "$HOME"`. Ako jedan aktivni račun traži više prostora ili RAM-a nego što poslužitelj daje, potreban je veći radni prostor ili drugi računalni resurs.

## Organizacija koda i izlaza

| Datoteka | Što u njoj potražiti |
|---|---|
| `postavke.json` | Šest temperatura, elektronska i fononska mreža, ravni valovi, pomaci i pragovi |
| `qe.py` | Geometrija monosloja, pretvorba temperature, tekst QE ulaza, izvođenje i provjera konvergencije |
| `studija.py` | Redoslijed relaksacije, DFPT-a i zamrznutih pomaka |
| `modovi.py` | Maseno ponderirani vektori, sektori savijanja i preklapanje modova |
| `izvjestaj.py` | Zaseban graf disperzije svake temperature i usporedni grafovi |

U svakoj temperaturnoj mapi `primjeri_ulaza` služi za čitanje prije izvođenja. **Stvarni ulazi** za fonone i pojedinačne amplitude nastaju uz završne geometrije u `fononi_1x1`, `gamma_2x2_ref`, `pes/<smjer>` i `gamma_2x2_rel`. U `relax_2x2` možeš usporediti `pocetna.extxyz` i `geometrija.extxyz`. U datoteci `gamma.npz` stupci matrice `e_massweighted` jesu svojstveni vektori, a niz `qe_e_complex` čuva izvorni QE izbor faze. Prethodni izračuni relaksacije, disperzije i Γ s istim postavkama mogu se ponovno koristiti. Ako imaš stari jedno-modni izlaz direktno u `pes/`, nove usporedbe čitaju imenovane podmape `pes/M1` i `pes/3Q`.

Istu vrijednost `degauss` možeš pronaći u svim elektronskim ulazima jedne temperature, uključujući relaksaciju i svaki zamrznuti pomak. `ph.x` preko `prefix` i `outdir` učitava elektronsko rješenje odgovarajućeg SCF-a i računa njegov odziv.

Kao u prvom dijelu projekta, `naprezanje_2d.csv` daje preostalo ravninsko naprezanje nakon ionske relaksacije. QE ispisuje tlak s pozitivnim predznakom za kompresiju. Za usporedbu s ASE-om kod promijeni predznak, pretvori Ry/bohr³ u eV/Å³ te pomnoži visinom ćelije. Ravninska konstanta ostaje zadanih 3,54 Å na svim temperaturama.

Početna mreža 24×24×1 služi za usporedbu svih temperatura. Za kvantitativan zaključak provjeri znak i veličinu najniže frekvencije na gušćim mrežama, posebice blizu eventualnog prijelaza. Za izmjenu napravi kopiju `postavke.json` i proslijedi je s `--config`. Različite konfiguracije dobivaju različite mape rezultata. Za test gušće fononske mreže možeš povećati `q_mesh` na 8 uz `k_mesh` 24 ili 48.

U osnovnom izvođenju neovisne grane dijele CPU redom. Na računalu s dovoljno memorije `python studija.py run --workers 2` može istodobno računati relaksiranu i simetričnu strukturu iste temperature. `TISE2_THREADS` postavlja broj OpenMP niti po procesu (zadano 2), a `TISE2_MPI_RANKS` i po potrebi `MPIEXEC` pokreću MPI inačice `pw.x` i `ph.x`. Valne funkcije koje više nisu potrebne brišu se iz pripadajućeg QE direktorija `tmp`. U slučaju pogreške ostaju ispis, ulaz i dijagnostika u `problemi.json`.

U prethodnom inicijalizacijskom ispisu (`nstep=0`) za stvarnu slučajno pomaknutu 2×2 pri 100 K, QE je pronašao 74 elektronske k-točke nakon redukcije i procijenio više od 6,40 GB dinamičke memorije po procesu. Zadani `--workers 1` obrađuje dvije skupe grane redom. Postavkom `--workers 2` mogu se izvršavati istodobno.
