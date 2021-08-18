from queue import Queue
import numpy as np
import ntpath
import os
from time import time

def combine_columns(singletonMatrixPartition):
    '''
    Funzione che combina due o più colonne nel calcolo del vettore rappresentativo di un insieme lamda
    composto da due o più elementi del dominio
    :param singletonMatrixPartition: partizione della matrice dei vettori rappresentativi dei vettori singoletti
    :return: il vettore rappresentativo dell'insieme lamda composto da due o più elementi di M
    '''
    #VIENE UTILIZZATO -1 al posto di 'x'
    res = np.zeros((singletonMatrixPartition.shape[0], 1), dtype=np.int64)
    while singletonMatrixPartition.shape[1] >= 2:
        c1 = np.array(singletonMatrixPartition[:, 0], dtype=np.int64)
        c2 = np.array(singletonMatrixPartition[:, 1], dtype=np.int64)
        i = 0
        while i < singletonMatrixPartition.shape[0]:
            if c1[i] != 0 and c2[i] != 0:
                res[i] = -1
            elif c1[i] == -1 or c2[i] == -1:
                res[i] = -1
            elif c1[i] != 0 and c2[i] == 0:
                res[i] = c1[i]
            elif c1[i] == 0 and c2[i] != 0:
                res[i] = c2[i]
            elif c1[i] == 0 and c2[i] == 0:
                res[i] = 0
            i += 1
        singletonMatrixPartition = np.delete(singletonMatrixPartition, [0, 1], 1)
        singletonMatrixPartition = np.concatenate((res, singletonMatrixPartition), axis=1)
    return singletonMatrixPartition

def build_representativeVector(lamda, singletonRepresentativeMatrix):
    '''
    Costruisce il vettore rappresentativo di un insieme lamda
    :param lamda: insieme di cui costruire il vettore rappresentative
    :param singletonRepresentativeMatrix: matrice dei vettori rappresentativi dei
    sottoinsiemi singoletti di M
    :return: il vettore rappresentativo
    '''
    if lamda.size == 0: #insieme vuoto
        return np.zeros((singletonRepresentativeMatrix.shape[0], lamda.size), dtype=np.int64)
    elif lamda.size == 1: #singoletto
        return np.array(singletonRepresentativeMatrix[:, lamda[0]-1])
    else: #insieme con almeno due elementi
        return combine_columns(singletonRepresentativeMatrix[:, lamda - 1])

def build_projection(lamda, representativeVector):
    '''
    Costruisce la proiezione del contenuto del
    vettore rappresentativo sull'insieme lamda associato
    :param lamda: insieme associato al vettore rappresentativo
    :param representativeVector: vettore rappresentativo di lamda
    :return: la proiezione intesa come insieme di elementi di lamda
    contenuti nel vettore rappresentativo
    '''
    projection = np.array([], dtype=np.int64)
    for elem in representativeVector:
        if elem in lamda:
            projection = np.append(projection,[elem])
    return set(projection)

def check(lamda, singletonRepresentativeMatrix):
    '''
    Effettua il controllo  CHECK sull'insieme lamda che gli viene passato
    :param lamda: insieme lamda di cui effettuare il CHECK
    :param singletonRepresentativeMatrix matrice dei vettori rappresentative dei singoletti
    :return: OK se l'insieme è un potenziale hitting set, KO se l'insieme non è un
    hitting set né può diventarlo, MHS se l'insieme è già un minimum hitting set
    '''

    representativeVector = build_representativeVector(lamda, singletonRepresentativeMatrix)
    projection = build_projection(lamda, representativeVector)

    #CHECK RULE
    if projection == set(lamda):
            if np.count_nonzero(representativeVector) == len(representativeVector):
                return 'MHS'
            else:
                return 'OK'
    else:
        return 'KO'

def output(lamda, counMHS, mapping):
    '''
    effettua l'ouput dell'insieme lamda, rivelatosi un mhs, insieme alla sua cardinalità e al conteggio corrente di mhs trovati
    :param lamda: di cui effettuare l'output
    :return: l'output di lamda mhs
    '''
    if mapping is None:
        print('MHS found: {} of dimension {}'.format(lamda, len(lamda)))
        print('MHS encountered : {} \n'.format(counMHS))
    else:
        lamda_remapped = [mapping[elem-1] for elem in lamda.tolist()]
        print('MHS found: {} of dimension {}'.format(lamda_remapped, len(lamda_remapped)))
        print('MHS encountered : {} \n'.format(counMHS))

def getMaps(indecesRemoved, MprimeLength):
    '''
    Funzione che calcola una mappa delle colonne di A dal dominio M al dominio M'
    quando si è effettuata una pre-elaborazione
    :param indecesRemoved: indici delle colonne rimosse nelle pre-elaborazione
    :param MprimeLength: lunghezza del nuovo dominio M'
    :return:
    '''
    indecesRemoved = np.array(indecesRemoved, dtype=np.int64)
    indecesRemoved = indecesRemoved + 1
    acc = 0
    mapping = []
    i = 1
    while len(mapping) < MprimeLength:
        while acc + i in indecesRemoved:
            acc += 1
        mapping.append(i+acc)
        i += 1
    return mapping

def getSingletonRepresentativeMatrix(A):
    '''
    Restituisce una matrice, intesa come vettore di vettori, contenente i vettori
    rappresentativi degli sotto-insiemi singoletti di M
    :param A: matrice di input
    :return: i vettori rappresentativi degli insiemi singoletti
    '''
    singletonMatrix = A
    for i in range(singletonMatrix.shape[0]):
        for j in range(singletonMatrix.shape[1]):
            if singletonMatrix[i, j] == 1:
                singletonMatrix[i, j] = j + 1
    return singletonMatrix

def mbase(A, timeEnabled=True, mapping=None):
    '''
    Algoritmo mbase per il calcolo di tutti i mhs data una matrice A che riassume il contenuto del dominio M
    e delle collezione di insieme N
    :param A: matrice di ingresso sopracitata
    :param timeEnabled: flag per il calcolo del tempo necessario alla computazione
    :param mapping: mapping delle colonne di A qualora si sia compiuta una pre-elaborazione
    '''
    if timeEnabled:
        start = time()
    coda = Queue(maxsize=0)
    coda.put(np.array([], dtype=np.int64))
    countMHS = 0

    singletonRepresentativeMatrix = getSingletonRepresentativeMatrix(A)
    M = np.array(list(range(1,A.shape[1]+1))) #prendo gli elementi di M
    #NB: M è già ordinato in ordine crescente per costruzione

    while not coda.empty():
        alpha = coda.get()

        if len(alpha) == 0:
            e = np.amin(M)
        else:
            e = np.amax(alpha) + 1
        while e <= np.amax(M):
            lamda = np.append(alpha, np.array(e, dtype=np.int64))
            #print('Esaminando lamda {}'.format(lamda))
            result = check(lamda, singletonRepresentativeMatrix)
            if result == 'OK' and e != np.amax(M):
                coda.put(lamda)
            elif result == 'MHS':
                countMHS += 1
                output(lamda, countMHS, mapping)
            # else:
            #     print('{} KO'.format(lamda))
            e += 1 #succ(e)
    if timeEnabled:
        end = time()
        print("MBASE required %.4f seconds to execute" % (end-start))

def getMatrixFromFile(filename):
    '''
    Legge un file .matrix dato il filename e restituisce una matrice a partire
    dal suo contenuto
    Per esempio a partire da './benchmarks1/74L85.000.matrix' si analizza il
    contenuto del file 74L85.000.matrix per creare la matrice
    :param filename: nome del file .matrix (indirizzo completo)
    :return: la matrice A costruita a partire dal contenuto del file
    '''
    print(ntpath.basename(filename))
    if os.stat(filename).st_size == 0:
        return np.array([[]],dtype=np.int64)
    else:
        file = open(filename, 'r')
        lines = file.readlines()
        lines = lines[5:]
        lines = [str.replace(line, ' -\n', '') for line in lines]
        lines = [str.replace(line, ' ', ',') for line in lines]
        #adesso abbiamo le righe della matrice A e possiamo costruirla
        rows = []
        for line in lines:
            rows.append(np.fromstring(line, sep=','))
        A = np.array(rows, dtype=np.int64)
        file.close()
        return A

def del_rows(A):
    '''
    Cancellazione di righe della matrice A che specificano per insiemi che sono superinsiemi di altre righe
    Per esempio, [1, 1, 1, 0] è superinsieme di [1, 1, 0, 0] perché la prima riga contiene tutti gli elementi che
    sono contenuti anche nella seconda (ne contiene perfino di più)
    Il ragionamento alla base del funzionamento è indicato tramite i seguenti esempi:
    Esempio 1:
    A 1 1 1 0 -
    B 0 1 1 0 =
      1 0 0 0  -> A super B
    Esempio 2:
    A 1 1 1 0 -
    B 0 1 1 1 =
      1 0 0 -1 -> A non super B
    Se c'è anche solo un -1 nel risultato significa che A non possiede almeno un elemento che invece B possiede. Ergo A non può essere superinsieme di B.
    Se invece la differenza consiste solo di '1' e '0' allora significa che A possiede tutti gli elementi che anche B possiede e forse anche elementi in più.
    Quindi si può pensare che se diff non contiene alcun '-1' allora A è superinsieme di B e quindi A va eliminato dalla matrice.
    Esempio 3:
    A 0 1 0 0
    B 0 0 1 0
      0 1 -1 0
    C'è un -1 quindi A non è superinsieme di B. Infatti si vede benissimo che A e B sono proprio disgiunti.
    Esempio 4:
    A 0 1 1 0
    B 1 1 1 0
     -1 0 0 0
    Contiene un -1 quindi A non è superinsieme di B ma vale il contrario, i.e. B è superinsieme di A
    Per come è fatto l'algoritmo B non verrebbe eliminato perché la ricerca si concentra sui casi in cui A è superinsieme di B.
    E' quindi necessario scorrere la matrice di ingresso dalla fine verso l'inizio così che si invertano i ruoli nella differenza vettoriale
    Sorge quindi un ulteriore problematica: in caso di righe duplicate (cioè A incluso in B e viceversa), l'indice delle righe da eliminare
    verrebbe indicato due volte. E' sufficiente quindi eliminare i duplicati dell'array degli indici delle righe da rimuovere prima di passare alla rimozione vera e propria
    :param A: matrice di input dell'algoritmo mhs
    :return: una matrice privata delle righe che sono superinsieme di altre irghe
    '''
    toBeRemoved = np.array([], dtype=np.int64)
    i, ii = 0, A.shape[0]-1
    while i <= A.shape[0] - 2 and ii >= 1:
        j, jj = i + 1, ii - 1
        while j <= A.shape[0] - 1 and jj >= 0:
            diff = A[i] - A[j]
            if np.count_nonzero(diff == -1) == 0:  # A[i] super A[j]
                toBeRemoved = np.append(toBeRemoved, [i])
            diff = A[ii] - A[jj]
            if np.count_nonzero(diff == -1) == 0:  # A[ii] super A[jj]
                toBeRemoved = np.append(toBeRemoved, [ii])
            j += 1
            jj -= 1
        i += 1
        ii -= 1
    toBeRemoved = np.unique(toBeRemoved)
    return np.delete(A, toBeRemoved, axis=0)

def del_cols(A):
    '''
    Cancellazione delle colonne interamente nulla
    Operazione da eseguirsi solamente dopo aver invocato del_cols
    Non viene eseguito alcun controllo sulla verifica di effettivo soddisfacimento di tale precedenza
    :param A:
    :return:
    '''
    zeroCols =(~np.all(A==0, axis=0)).tolist()
    indecesRemoved = []
    for i in range(len(zeroCols)):
        if zeroCols[i] == False:
            indecesRemoved.append(i)
    print('Columns dropped in preprocessing: ', [idx+1 for idx in indecesRemoved])
    return (A[:, zeroCols], indecesRemoved)

def pre_processing(A):
    '''
    Restituisce una tupla contentene la matrice A di N' righe e M' colonne a seguito della
    pre-elaborazione e gli indici delle colonne che sono state rimosse
    :param A: matrice A input dell'algoritmo
    '''
    start = time()
    Aprime = del_cols(del_rows(A))
    end = time()
    print("Preprocessing required %.6f seconds to execute" % (end - start))
    return Aprime

A = getMatrixFromFile(filename='c432.000.matrix')
if np.size(A) == 0:
    print("The specified file is empty. Can't start the computation")
else:
    (A_post, indecesRemoved) = pre_processing(A)
    #print(A_post, '\n')
    print('(|N| = {}, |M| = {}) \n'.format(A_post.shape[0], A_post.shape[1]))
    mbase(A_post, mapping=getMaps(indecesRemoved, A_post.shape[1]))
    #secondo giro con mbase normale
    print('(|N| = {}, |M| = {}) \n'.format(A.shape[0], A.shape[1]))
    mbase(A)

#ESEMPIO VISTO IN AULA -> OK!
    # A = np.matrix([[1, 1, 1, 0, 0, 0],
    #               [0, 1, 1, 1, 0, 1],
    #               [0, 1, 1, 0, 0, 0],
    #               [1, 0, 1, 0, 0, 1],
    #               [1, 1, 1, 1, 0, 0],
    #               [0, 1, 1, 1, 0, 1],
    #               [0, 0, 1, 1, 0, 0]], dtype=np.int64)
