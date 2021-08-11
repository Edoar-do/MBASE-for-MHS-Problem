from queue import Queue
import numpy as np
def build_representativeVector(lamda, singletonRepresentativeMatrix):
    '''
    Costruisce il vettore rappresentativo di un insieme lamda
    :param lamda: insieme di cui costruire il vettore rappresentative
    :param singletonRepresentativeMatrix: matrice dei vettori rappresentativi dei
    sottoinsiemi singoletti di M
    :return: il vettore rappresentativo
    '''
    #TODO

def build_projection(lamda, representativeVector):
    '''
    Costruisce la proiezione del contenuto del
    vettore rappresentativo sull'insieme lamda associato
    :param lamda: insieme associato al vettore rappresentativo
    :param representativeVector: vettore rappresentativo di lamda
    :return: la proiezione intesa come insieme di elementi di elementi di lamda
    contenuti nel vettore rappresentativo
    '''
    #TODO

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
    if set(projection) == set(lamda):
            if np.count_nonzero(representativeVector) == len(representativeVector):
                return 'MHS'
            else:
                return 'OK'
    else:
        return 'KO'


def output(lamda, dim):
    '''
    effettua l'ouput dell'insieme lamda, rivelatosi un mhs
    :param lamda: di cui effettuare l'output
    :param dim: dimensione della collezione N che è anche la dimensione del
    vettore rappresentativo di lamda
    :return: l'output di lamda mhs
    '''
    #vedere se stampare e basta oppure se creare una lista
    #che verrà restituita all'utente. In qualche modo devo stampare il numero di MHS
    #trovati e la cardinalità minima e massima di tali MHS

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
        lamda = coda.get()

        if len(lamda) == 0:
            e = min(M)
        else:
            e = max(lamda) + 1

        while e <= max(M):
            lamda.extend([e])
            result = check(lamda, singletonRepresentativeMatrix)
            if result == 'OK' and e != max(M):
                coda.put(lamda)
            elif result == 'MHS':
                output(lamda)
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
    A = np.matrix(rows, dtype=np.bool) #metto come tipo bool per risparmiare spazio
    file.close()
    return A


