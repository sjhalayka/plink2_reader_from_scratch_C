#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>



// Structure to represent the Plink2Reader object
typedef struct {
    FILE* pgen_file;
    FILE* pvar_file;
    FILE* psam_file;
    uint32_t variant_count;
    uint32_t sample_count;
    uint64_t file_size;
} Plink2Reader;

// Function to initialize and open files
void Plink2Reader_open(Plink2Reader* reader, const char* pgen_path, const char* pvar_path, const char* psam_path) {
    reader->pgen_file = fopen(pgen_path, "rb");
    reader->pvar_file = fopen(pvar_path, "r");
    reader->psam_file = fopen(psam_path, "r");

    if (!reader->pgen_file || !reader->pvar_file || !reader->psam_file) {
        fprintf(stderr, "Failed to open one or more PLINK2 files\n");
        exit(EXIT_FAILURE);
    }
}

// Function to close all files
void Plink2Reader_close(Plink2Reader* reader) {
    if (reader->pgen_file) fclose(reader->pgen_file);
    if (reader->pvar_file) fclose(reader->pvar_file);
    if (reader->psam_file) fclose(reader->psam_file);
}

// Function to read the header and metadata
void Plink2Reader_readHeader(Plink2Reader* reader) {
    char magic[2];
    fread(magic, 1, 2, reader->pgen_file);

    if (magic[0] != 0x6c || magic[1] != 0x1b) {
        fprintf(stderr, "Invalid PGEN file format\n");
        exit(EXIT_FAILURE);
    }

    char storage_mode;
    fread(&storage_mode, 1, 1, reader->pgen_file);

    if (storage_mode != 0x10) {
        fprintf(stderr, "Unsupported storage mode\n");
        exit(EXIT_FAILURE);
    }

    fread(&reader->variant_count, sizeof(uint32_t), 1, reader->pgen_file);
    fread(&reader->sample_count, sizeof(uint32_t), 1, reader->pgen_file);

    fseek(reader->pgen_file, 0, SEEK_END);
    reader->file_size = ftell(reader->pgen_file);
    fseek(reader->pgen_file, 11, SEEK_SET); // Set to start of genotype data
}

// Function to read a chunk of genotypes
void Plink2Reader_readGenotypesChunk(Plink2Reader* reader, int** genotypes, uint32_t start_variant, uint32_t end_variant, uint32_t start_sample, uint32_t end_sample) {
    if (end_variant >= reader->variant_count || end_sample >= reader->sample_count) {
        fprintf(stderr, "Requested chunk is out of range\n");
        exit(EXIT_FAILURE);
    }

    uint32_t num_variants = end_variant - start_variant;
    uint32_t num_samples = end_sample - start_sample;

    genotypes = (int**)malloc(num_samples * sizeof(int*));
    for (uint32_t i = 0; i < num_samples; i++) {
        genotypes[i] = (int*)malloc(num_variants * sizeof(int));
    }

    uint64_t start_pos = 11 + (start_variant * reader->sample_count + start_sample) / 4;
    fseek(reader->pgen_file, start_pos, SEEK_SET);

    uint32_t bytes_to_read = (end_variant - start_variant) * (end_sample - start_sample);
    uint8_t* file_chunk = (uint8_t*)malloc(bytes_to_read);

    fread(file_chunk, 1, bytes_to_read, reader->pgen_file);

    uint32_t file_chunk_index = 0;

    for (uint32_t variant = start_variant; variant < end_variant; variant++) {
        for (uint32_t sample = start_sample; sample < end_sample; sample++) {
            int genotype = file_chunk[file_chunk_index] & 0x03;
            genotypes[sample - start_sample][variant - start_variant] = (genotype == 3) ? -1 : genotype;
            file_chunk_index++;
        }
    }

    free(file_chunk);
}

// Function to read a chunk of variant info
void Plink2Reader_readVariantInfoChunk(Plink2Reader* reader, char** variant_ids, uint32_t start_variant, uint32_t end_variant) {
    if (end_variant >= reader->variant_count) {
        fprintf(stderr, "Requested chunk is out of range\n");
        exit(EXIT_FAILURE);
    }

    char line[1024];
    fgets(line, sizeof(line), reader->pvar_file); // Skip header line

    for (uint32_t i = 0; i < start_variant; i++) {
        fgets(line, sizeof(line), reader->pvar_file); // Skip to the start variant
    }

    for (uint32_t i = start_variant; i < end_variant; i++) {
        fgets(line, sizeof(line), reader->pvar_file);
        variant_ids[i - start_variant] = strdup(strtok(line, "\t"));
    }
}

// Function to read a chunk of sample info
void Plink2Reader_readSampleInfoChunk(Plink2Reader* reader, char** sample_ids, uint32_t start_sample, uint32_t end_sample) {
    if (end_sample >= reader->sample_count) {
        fprintf(stderr, "Requested chunk is out of range\n");
        exit(EXIT_FAILURE);
    }

    char line[1024];
    fgets(line, sizeof(line), reader->psam_file); // Skip header line

    for (uint32_t i = 0; i < start_sample; i++) {
        fgets(line, sizeof(line), reader->psam_file); // Skip to the start sample
    }

    for (uint32_t i = start_sample; i < end_sample; i++) {
        fgets(line, sizeof(line), reader->psam_file);
        sample_ids[i - start_sample] = strdup(strtok(line, "\t"));
    }
}

int main(void) {
    Plink2Reader reader;

    Plink2Reader_open(&reader, "plink2.pgen", "plink2.pvar", "plink2.psam");
    Plink2Reader_readHeader(&reader);

    printf("Variant count: %u\n", reader.variant_count);
    printf("Sample count: %u\n", reader.sample_count);

    Plink2Reader_close(&reader);
    return 0;
}
