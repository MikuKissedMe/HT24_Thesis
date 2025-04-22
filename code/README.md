Du bör konstruera några enkla egna exempel. Som exempel, hur visualiseras en kolumn med identiska
bokstäver? Dvs en kolumn som AAAAAAA och vilket träd som helst? Och jämför med en kolumn där varje
position skiljer sig åt, som ”ARVESTID” :-)

Ett annat exempel är kolumnen AAAACCCC, och med träd av olika slag, tex det balanserade trädet
(((a,b),(c,d)), ((e,f), (g,h))) och jämfört med det helt obalanserade trädet (a, (b, (c, (d, (e, (f,
(g, h))))))).

Prova även ACCCCCC med de två träden. Här borde ju PIC göra de obalanserade trädet ge en helt annan
visualisering än det balanserade, eller hur?

# Data

## ex1

* Två första kolumnerna är helt konserverade: bara "A".
* Sen kommer två kolumner med fullständig spridning av "ARVESTID".
* Kolumnerna 5 och 6 är Hälften vardera av A och C på övre respektive undre halvan.
* Kolumnerna 7 och 8 har omväxlande A och C på varannan rad.
* Kolumnerna 9 och 10 har "R" överallt förutom första (a) respektive sista (h) sekvensen.
  - I detta fall borde utfallet bli mycket olika beroende på vilket träd man använder.

### Träd

* ex1_t1.tree: Helt balanserat träd.
* ex1_t2.tree: Helt obalanserat, med "a" som eget delträd.
* ex1_t3.tree: Helt obalanserat, som ex1_t2, men med "h" istf "a" som eget delträd.



$ for x in 1 2 3; do nw_display ex1_t$x".tree"; done
                                                    +------------------------+ a
                          +-------------------------+
                          |                         +------------------------+ b
 +------------------------+
 |                        |                         +------------------------+ c
 |                        +-------------------------+
 |                                                  +------------------------+ d
=|
 |                                                  +------------------------+ e
 |                        +-------------------------+
 |                        |                         +------------------------+ f
 +------------------------+
                          |                         +------------------------+ g
                          +-------------------------+
                                                    +------------------------+ h

 |------------|-----------|------------|------------|-----------|-------------
 0         0.05         0.1         0.15          0.2        0.25
 substitutions/site

 +---------------------------------------------------------------------------+ a
 |
=|          +----------------------------------------------------------------+ b
 |          |
 +----------+          +-----------------------------------------------------+ c
            |          |
            +----------+          +------------------------------------------+ d
                       |          |
                       +----------+         +--------------------------------+ e
                                  |         |
                                  +---------+          +---------------------+ f
                                            |          |
                                            +----------+          +----------+ g
                                                       +----------+
                                                                  +----------+ h

 |---------------------|--------------------|---------------------|-----------
 0                   0.1                  0.2                   0.3
 substitutions/site

 +---------------------------------------------------------------------------+ h
 |
=|          +----------------------------------------------------------------+ b
 |          |
 +----------+          +-----------------------------------------------------+ c
            |          |
            +----------+          +------------------------------------------+ d
                       |          |
                       +----------+         +--------------------------------+ e
                                  |         |
                                  +---------+          +---------------------+ f
                                            |          |
                                            +----------+          +----------+ g
                                                       +----------+
                                                                  +----------+ a

 |---------------------|--------------------|---------------------|-----------
 0                   0.1                  0.2                   0.3
