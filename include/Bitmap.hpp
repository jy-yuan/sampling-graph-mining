#ifndef BITMAP_HPP
#define BITMAP_HPP

#define WORD_OFFSET(i) ((i) >> 6)
#define BIT_OFFSET(i) ((i) & 0x3f)

#include <stddef.h>

class Bitmap {
	private:
		size_t size;
		unsigned long * data;
	public:
		Bitmap();
		Bitmap(size_t size_); // allocate an empty bitmap
		~Bitmap();

		void clear();
		void fill();
		unsigned long get_bit(size_t i);
		void set_bit(size_t i);
};

#endif
