# Sadržaj DFT paketa

Vodič i skripte ažurirani su 23. studenoga 2026. Za početak otvori **Vodic_DFT_TiSe2.pdf**. Naredbe za jednu temperaturu i cijelu seriju nalaze se na zadnjoj stranici PDF-a i u datoteci **README.md**.

Priložena mapa `racuni` čuva prethodnu snimku od 22. studenoga 2026. u 16:32 po zagrebačkom vremenu. Osnovne postavke imaju oznaku `k24_q6_17503ad3bf`, a tehnička provjera pri 1600 K na gruboj mreži oznaku `k4_q2_b4968b2f3c`. U svakoj je mapi sačuvan tadašnji `rezultati/PREGLED.md`. Ta snimka prethodi uvođenju usporedbe M₁ i 3Q pa njezine stare tablice stanja imaju jedan stupac PES.

Nakon tvog pokretanja `studija.py run` skripta za svaku dovršenu temperaturu sprema vektore M₁, M₂, M₃ i njihovu jednako ponderiranu kombinaciju te dvije datoteke `potencijal.csv`, jednu u `pes/M1`, drugu u `pes/3Q`. `studija.py report` osvježava tablice iz već dovršenih računa. Priloženi brojevi za minimum oba presjeka u vodiču očitani su iz MACE dijela bilježnice. Novi DFT presjeci računaju se pri tvom pokretanju.
