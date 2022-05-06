#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int cmpfunc(const void* a, const void* b);
int readInput(char *filename, int **array);
int writeOutput(char *fileName, int *output, int arraySize);
int checkResult(int* result, int arraySize);
int findMedian(int* array, int arraySize);
void merge(int* A, int* B, int arraySizeA, int arraySizeB, int** result);
int strategyIsValid(int strategy);
int meanOfArray(int* array, int arraySize);

int main(int argc, char* argv[]) {
    char *inputName = argv[1];
    char *outputName = argv[2];
    int chosenStrategy = atoi(argv[3]);

    int size, cubeRank;
    int worldRank;
    int groupRank;
    int groupSize;
    int blockSize;

    int* input;
    int* output;

    int arraySize;

    int pivot;

    MPI_Status status;
    MPI_Comm nCube, groupComm;
    MPI_Request send_request_A, recv_request_A, send_request_B, recv_request_B;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    int elementCountPerProcess[size];
    int currentDisplacement = 0;
    int scatterDisplacement[size];
    int finalElementCount[size];
    int displacementArray[size];

    if (worldRank == 0)
    {
        // Check if chosen strategy is valid
        if (!strategyIsValid(chosenStrategy)) {
            return 2;
        }

        if (0 > (arraySize = readInput(inputName, &input)))
            return 2;

        // If arraySize is divisible by size, place an equal amount of elements in each process' array
        if (arraySize%size == 0) {
            blockSize = arraySize / size;

            for (int i = 0; i < size; i++) {

                elementCountPerProcess[i] = blockSize;
                scatterDisplacement[i] = currentDisplacement;
                currentDisplacement += blockSize;

            }
        }

        // Place an extra element in the arraySize%size first process arrays, then arraySize/size in the rest
        else {
            int elementsLeft = arraySize;
            int elementsPerProcess = arraySize / size;

            for (int i = 0; i < size; i++) {
                if (i < (arraySize%size)) {
                    elementCountPerProcess[i] = elementsPerProcess + 1;
                    scatterDisplacement[i] = currentDisplacement;
                    currentDisplacement += elementsPerProcess+1;
                    elementsLeft -= elementsPerProcess-1;
                }
                else {
                    elementCountPerProcess[i] = elementsPerProcess;
                    scatterDisplacement[i] = currentDisplacement;
                    currentDisplacement += elementsPerProcess;
                    elementsLeft -= elementsPerProcess;
                }
            }
        }

        output = malloc(sizeof(int)*arraySize);
    }

    MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Create hypercube
    int dimensionAmount = (int) log2(size);
    int dims[dimensionAmount];
    int period[dimensionAmount];
    int coords[dimensionAmount];

    for (int i = 0; i < dimensionAmount; i++) {
        period[i] = 1;
        dims[i] = 2;
    }

    MPI_Cart_create(MPI_COMM_WORLD, dimensionAmount, dims, period, 1, &nCube);
    MPI_Comm_dup(nCube, &groupComm);
    MPI_Comm_size(groupComm, &groupSize);
    MPI_Comm_rank(nCube, &cubeRank);
    MPI_Cart_coords(nCube, cubeRank, dimensionAmount, coords);

    // Broadcast the amounts of elements per process and allocate an appropriate amount of memory
    MPI_Bcast(elementCountPerProcess, size, MPI_INT, 0, MPI_COMM_WORLD);
    blockSize = elementCountPerProcess[cubeRank];
    int* localBlock = malloc(blockSize*sizeof(int));

    MPI_Scatterv(input, elementCountPerProcess, scatterDisplacement, MPI_INT, localBlock, blockSize, MPI_INT, 0, nCube);

    int** sendArray = malloc(sizeof(int*) * 2);
    sendArray[0] = malloc(sizeof(int) * arraySize);
    sendArray[1] = malloc(sizeof(int) * arraySize);
    int* recvArray = malloc(arraySize * sizeof(int));

    int* sendAmount = malloc(2*sizeof(int));

    // Start the timer
    double timeStart = MPI_Wtime();

    // Perform initial local sort
    qsort(localBlock, (size_t)blockSize, sizeof(int), cmpfunc);

    for (int d = 0; d < dimensionAmount; d++) {

        // Calculate new median and broadcast it to group
        switch (chosenStrategy)
        {
            case 1:
            {
                if (groupRank == 0) {
                    pivot = findMedian(localBlock, blockSize);
                }
                break;
            }
            case 2:
            {
                // Find median of each process' array
                pivot = findMedian(localBlock, blockSize);
                int* medianArray = malloc(sizeof(int) * groupSize);
                
                MPI_Gather(&pivot, 1, MPI_INT, medianArray, 1, MPI_INT, 0, groupComm);

                // Determine the median of all medians
                if (groupRank == 0) {

                    qsort(medianArray, (size_t)groupSize, sizeof(int), cmpfunc);
                    pivot = findMedian(medianArray, groupSize);

                }
                free(medianArray);
            }

            case 3:
            {
                // Find median of each process' array
                pivot = findMedian(localBlock, blockSize);
                int* medianArray = malloc(sizeof(int) * groupSize);
                
                MPI_Gather(&pivot, 1, MPI_INT, medianArray, 1, MPI_INT, 0, groupComm);
                
                // Find mean of all medians
                if (groupRank == 0) {
                    pivot = meanOfArray(medianArray, groupSize);
                }
                free(medianArray);
            }
        }

        MPI_Bcast(&pivot, 1, MPI_INT, 0, groupComm);

        // Divide local array into upper and lower parts
        sendAmount[0] = 0; //Upper
        sendAmount[1] = 0; //Lower
        for (int i = 0; i < blockSize; i++) {
            if (localBlock[i] > pivot) {
                sendArray[0][sendAmount[0]] = localBlock[i];
                sendAmount[0] += 1;
            }
            else {
                sendArray[1][sendAmount[1]] = localBlock[i];
                sendAmount[1] += 1;
            }
        }

        MPI_Comm_split(groupComm, coords[d], groupRank, &groupComm);
        MPI_Comm_size(groupComm, &groupSize);
        MPI_Comm_rank(groupComm, &groupRank);

        // Determine group rank and nearest neighbor
        int source, destination;
        int receiveSize;

        // If senderCoordinate is 0, then receiverCoordinate is 1.
        int senderCoordinate = coords[d];
        int receiverCoordinate = (coords[d]+1)%2;

        MPI_Cart_shift(nCube, d, 1, &source, &destination);

        // Send the upper or lower part of the local array to destination
        MPI_Sendrecv(&(sendAmount[senderCoordinate]), 1, MPI_INT, destination, 0, &receiveSize, 1, MPI_INT, destination, 0, nCube, &status);
        MPI_Sendrecv(&(sendArray[senderCoordinate][0]), sendAmount[senderCoordinate], MPI_INT, destination, 1, &(recvArray[0]), receiveSize, MPI_INT, destination, 1, nCube, &status);
        
        // Allocate more memory if necessary, or reduce size of local array
        if (blockSize < sendAmount[receiverCoordinate] + receiveSize) {
            localBlock = realloc(localBlock, (sendAmount[receiverCoordinate]+receiveSize)*sizeof(int));
            blockSize = sendAmount[receiverCoordinate]+receiveSize;
        }
        else if (blockSize > (sendAmount[receiverCoordinate] + receiveSize)) blockSize = sendAmount[receiverCoordinate]+receiveSize;

        // Merge the received array with the array the process kept 
        merge(sendArray[receiverCoordinate], recvArray, sendAmount[receiverCoordinate], receiveSize, &localBlock);

    }

    double maxTime;
    double executionTime = MPI_Wtime() - timeStart;

    // Find the amount of elements per process
    MPI_Gather(&blockSize, 1, MPI_INT, finalElementCount, 1, MPI_INT, 0, nCube);

    displacementArray[0] = 0;
    int displacement = finalElementCount[0];

    if (cubeRank == 0) {
        for (int i = 1; i < size; i++) {
            displacementArray[i] = displacement;
            displacement += finalElementCount[i];
        }
    }
    MPI_Gatherv(localBlock, blockSize, MPI_INT, output, finalElementCount, displacementArray, MPI_INT, 0, nCube);
    // Find the largest execution time
    MPI_Reduce(&executionTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (worldRank == 0) {
        writeOutput(outputName, output, arraySize);
        printf("%lf\n", maxTime);
        free(input);
        free(output);
    }

    free(localBlock);
    free(sendArray[0]);
    free(sendArray[1]);
    free(sendArray);
    free(sendAmount);
    free(recvArray);
    MPI_Comm_free(&nCube);
    MPI_Comm_free(&groupComm);
    MPI_Finalize();
}

int meanOfArray(int* array, int arraySize) {
    int sum = 0;
    for (int i = 0; i < arraySize; i++) {
        sum += array[i];
    }
    return sum / arraySize;
}

int strategyIsValid(int strategy) {
    switch (strategy)
    {
        case 1:
            return 1;
        case 2:
            return 1;
        case 3:
            return 1;
        default:
            printf("Invalid strategy.\n");
            return 0;
    }
}

// Helper function for qsort
int cmpfunc(const void* a, const void* b) {
    return ( *(int*)a - *(int*)b );
}

int findMedian(int* array, int arraySize) {
    int median = 0;
    if (arraySize%2 == 0) {
        median = (array[(arraySize-1)/2] + array[arraySize/2])/2;
    }

    else median = array[arraySize/2];

    return median;
}

// Function used to merge two sorted lists
void merge(int* A, int* B, int arraySizeA, int arraySizeB, int** result) {
    int i = 0, j = 0, k = 0;

    while (i < arraySizeA && j < arraySizeB) {
        if (A[i] < B[j]) {
            (*result)[k] = A[i];
            k++; i++;
        }
        else {
            (*result)[k] = B[j];
            k++; j++;
        }
    }

    while (i < arraySizeA) {
        (*result)[k] = A[i];
        k++; i++;
    }

    while (j < arraySizeB) {
        (*result)[k] = B[j];
        k++; j++;
    }
}

// Function for handling input. Inspired by the readInput function given in A1.
int readInput(char *filename, int **array)
{
    FILE *file;
    if (NULL == (file = fopen(filename, "r")))
    {
        printf("Error opening file.\n");
        return -1;
    }
    int arraySize;
    if (EOF == fscanf(file, "%d", &arraySize))
    {
        perror("Couldn't read matrix size from input file");
        return -1;
    }

    *array = malloc(sizeof(int)*arraySize);

    for (int j = 0; j < arraySize; j++)
    {
        if (EOF == fscanf(file, "%d", &((*array)[j])))
        {
            perror("Couldn't read elements from input file");
            return -1;
        }
    }

    if (0 != fclose(file))
    {
        perror("Warning: couldn't close input file");
    }
    return arraySize;
}

// Function for handling output. Inspired by the readInput function given in A1.
int writeOutput(char *fileName, int *output, int arraySize) {
    FILE *file;
    if (NULL == (file = fopen(fileName, "w"))) {
        perror("Couldn't open output file");
        return -1;
    }
    for (int j = 0; j < arraySize; j++) {
        if (0 > fprintf(file, "%d ", output[j])) {
            perror("Couldn't write to output file");
        }
    }
    if (0 > fprintf(file, "\n")) {
        perror("Couldn't write to output file");
    }
    if (0 != fclose(file)) {
        perror("Warning: couldn't close output file");
    }
    return 0;
}

/*
int checkResult(int* result, int arraySize) {
    for (int i = 1; i < arraySize; i++) {
        if (result[i] < result[i-1]) return 0;
    }
    return 1;
}
*/