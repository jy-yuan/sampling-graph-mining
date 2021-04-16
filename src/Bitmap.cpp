#include "Bitmap.hpp"

Bitmap::Bitmap(): size(0), data(NULL) {
}
		
Bitmap::Bitmap(size_t size_): size(size_) {
	data = new unsigned long [WORD_OFFSET(size) + 1];
	clear();
}
		
Bitmap::~Bitmap() {
	delete [] data;
}

void Bitmap::clear() {
	size_t bm_size = WORD_OFFSET(size);
#pragma omp parallel for 
	for (size_t i = 0; i <= bm_size; ++ i) {
		data[i] = 0;
	}
}

void Bitmap::fill() {
	size_t bm_size = WORD_OFFSET(size);
#pragma omp parallel for 
	for (size_t i = 0; i <= bm_size; ++ i) {
		data[i] = 0xffffffffffffffff;
	}
	data[bm_size] = 0;
	for (size_t i = (bm_size << 6); i < size; ++ i) {
		data[bm_size] |= 1ul << BIT_OFFSET(i);
	}
}

unsigned long Bitmap::get_bit(size_t i) {
	return data[WORD_OFFSET(i)] & (1ul << BIT_OFFSET(i));
}

void Bitmap::set_bit(size_t i) {
	__sync_fetch_and_or(data + WORD_OFFSET(i), 1ul << BIT_OFFSET(i));
}

