from queue import Queue
import numpy as np

def combine_columns(singletonMatrixPartition):
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
    if len(lamda) == 0:
        return np.zeros((singletonRepresentativeMatrix.shape[0], len(lamda)), dtype=np.int64)
    elif len(lamda) == 1: #singoletto
        return np.array(singletonRepresentativeMatrix[:, lamda[0]+1])
    else:
        temp = np.zeros((singletonRepresentativeMatrix.shape[0], len(lamda)), dtype=np.int64)
        for i in range(len(lamda)):
            temp[:, i] = singletonRepresentativeMatrix[:, lamda[i]+1]
        return combine_columns(temp)

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


def output(lamda):
    '''
    effettua l'ouput dell'insieme lamda, rivelatosi un mhs
    :param lamda: di cui effettuare l'output
    :return: l'output di lamda mhs
    '''
    #vedere se stampare e basta oppure se creare una lista
    #che verrà restituita all'utente. In qualche modo devo stampare il numero di MHS
    #trovati e la cardinalità minima e massima di tali MHS

    #VERSIONE DI PROVA
    print('MHS trovato: {} di cardinalità {}'.format(lamda, len(lamda)))
    print()

def getSingletonRepresentativeMatrix(A):
    '''
    Restituisce una matrice, intesa come vettore di vettori, contenente i vettori
    rappresentativi degli sotto-insiemi singoletti di M
    :param A: matrice di input
    :return: i vettori rappresentativi degli insiemi singoletti
    '''
    singletonMatrix = A
    i = 0
    while i < singletonMatrix.shape[0]:
        j = 0
        while j < singletonMatrix.shape[1]:
            if singletonMatrix[i,j] == 1:
                singletonMatrix[i,j] = j+1
            j += 1
        i += 1
    return singletonMatrix

def mbase(A):
    coda = Queue(maxsize=0)
    coda.put([])

    singletonRepresentativeMatrix = getSingletonRepresentativeMatrix(A)
    M = list(range(1,A.shape[1]+1)) #prendo gli elementi di M
    #NB: M è già ordinato in ordine crescente per costruzione

    while not coda.empty():
        alpha = coda.get()

        if len(alpha) == 0:
            e = min(M)
        else:
            e = max(alpha) + 1
        while e <= max(M):
            lamda = alpha + [e]
            print('Esaminando lamda {}'.format(lamda))
            result = check(lamda, singletonRepresentativeMatrix)
            if result == 'OK' and e != max(M):
                coda.put(lamda)
            elif result == 'MHS':
                output(lamda)
            else:
                print('{} KO'.format(lamda))
            e += 1 #succ(e)

def getMatrixFromFile(filename):
    '''
    Legge un file .matrix dato il filename e restituisce una matrice a partire
    dal suo contenuto
    Per esempio a partire da './benchmarks1/74L85.000.matrix' si analizza il
    contenuto del file 74L85.000.matrix per creare la matrice
    :param filename: nome del file .matrix (indirizzo completo)
    :return: la matrice A costruita a partire dal contenuto del file
    '''

    file = open(filename, 'r')
    lines = file.readlines()
    lines = lines[5:]
    lines = [str.replace(line, ' -\n', '') for line in lines]
    lines = [str.replace(line, ' ', ',') for line in lines] #trim
    #adesso abbiamo le righe della matrice A e possiamo costruirla
    rows = []
    for line in lines:
        rows.append(np.fromstring(line, sep=','))
    A = np.matrix(rows, dtype=np.int64)
    file.close()
    return A


A = getMatrixFromFile(filename='74L85.000.matrix')
print(A)
print()
mbase(A)