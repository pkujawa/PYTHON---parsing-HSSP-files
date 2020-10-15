import re
import math
# Otwieranie pliku o nazwie plik.txt
plik = open('plik.txt', 'r')
#text to cala zawartosc tekstowa pliku
text = plik.read()
# 2 zmienne odczytane z pierwszych linii w pliku, 'SEQLENGTH' oraz 'NALIGN' to dlugosc sekwencji oraz liczba sekwencji
seqlength = int([line for line in text.split("\n") if "SEQLENGTH" in line][0].split()[1].strip())
nalign = int([line for line in text.split("\n") if "NALIGN" in line][0].split()[1].strip())


def parse_hsp(text):
    """
    fukcja dzielaca plik na sekcje PROTEINS oraz ALIGNMENTS, sekcja PROTEINS konczy sie na pierwszej sekcji ALIGNMENTS,
    ostatnia sekcja ALIGNMENTS konczy sie na sekcji SEQUENCE
    return: funckja zwraca indeksy poczatku sekcji proteins, alignments oraz sequence
    """
    indeksyA = []
    proteins = re.search(r'^## PROTEINS', text, flags=re.MULTILINE)
    indeksP = proteins.start()
    alignments = re.finditer(r'^## ALIGNMENTS', text, flags=re.MULTILINE)
    # dla kazdej sekwencji ALIGMENTS znajdz poczatek sekcji
    for a in alignments:
        indeksyA.append(a.start())
    koniec = re.search(r'^## SEQUENCE', text, flags=re.MULTILINE)
    indeksK = koniec.start()
    return (indeksP, indeksyA, indeksK)


def parse_proteins(indeksy):
    """
    funkcja pobiera ID protein zapisne w sekcji PROTEINS i zapisuje je do tablicy IDs
    :param indeksy: indeksy z fukcji parse_hsp
    :return: tablic ID protein
    """
    IDs = []
    # sekcja proteins zawarta jest w tekscie od inceksy indeksP z parse_hsp do poczatku pierwszej sekcji alignments
    proteins = text[indeksy[0]:indeksy[1][0]]
    lines = proteins.split('\n')
    lines = lines[2:-1] #omijamy 2 pierwsze linie
    for line in lines:
        ID = line.split()[2] #id znajduje sie w 3 kolumnie
        IDs.append(ID)
    return IDs


def parse_alignments(indeksy):
    """
    funckja pobierjaca ze wszystkich sekcji ALIGNMENTS te fragmenty linii, ktore zawieraja sekwencje aminokwasowe
    :param indeksy: indeksy z fukcji parse_hsp
    :return:
    """
    sequences = []
    for x in range(len(indeksy[1]) - 1): #dla kazdej seksji ALIGNMENTS
        alignments = text[indeksy[1][x]:indeksy[1][x + 1]] #kolejna sekcja ALIGMENTS
        lines = alignments.split('\n')
        lines = lines[2:-1] #pomijane sa 2 pierwsze linie
        for line in lines:
            seq = line[51:121] #sekwencja zapisana jest od znaku 50 do 120
            seq = seq.replace(' ', '-') # przerwy zamieniane sa na '-'
            sequences.append(seq)

    # osobne postepowanie dla ostatniej sekcji ALIGMENTS, ktora konczy sie tam, gdzie zaczyna sie sekcja SEQUENCE
    last = text[indeksy[1][-1]:indeksy[2]]
    lines = last.split('\n')
    lines = lines[2:-1]
    for line in lines:
        seq = line[51:121]
        seq = seq.replace(' ', '-')
        sequences.append(seq)
    return sequences


def write_result(IDs, sequences):
    """

    :param IDs: lista ID protein
    :param sequences: lista linii z sekwencjami z sekcji ALIGNMENTS
    :return:
    """
    result = []
    times = math.ceil(nalign / 70) # obliczenie ile razy wystepuje sekcja ALIGNMENTS
    for x in range(times): # dla kazdej sekcji ALIGNMENTS
        min = x * seqlength #poczatek sekcji
        max = (x + 1) * seqlength #koniec sekcji
        for letter in range(0, 70):
            s = ''
            # scal n-ty znak z 70 kolejnych linii
            for element in range(min, max):
                s += str(sequences[element][letter])
            result.append(s)
    # wypisz dane w formacie 'ID: sekwencja'
    for ID in range(0, len(IDs)):
        print (IDs[ID] + ": " + result[ID])


indeksy = parse_hsp(text)
write_result(parse_proteins(indeksy), parse_alignments(indeksy))
