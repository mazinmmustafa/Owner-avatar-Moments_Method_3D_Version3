//
#include "file.hpp"

file_t::file_t(){

}

file_t::~file_t(){
}

void file_t::open(const char *filename, const char mode){
    assert_error(this->is_open==false, "file is already open");
    this->filename = (char*)calloc(this->max_length, sizeof(char));
    assert(this->filename!=null);
    strcpy(this->filename, filename);
    assert_error(mode=='r'||mode=='w'||mode=='a', "invalid file mode");
    this->mode = mode;
    switch (this->mode){
    case 'w':
        this->file_ptr = fopen(this->filename, "w");
        break;
    case 'r':
        this->file_ptr = fopen(this->filename, "r");
        break;
    case 'a':
        this->file_ptr = fopen(this->filename, "a");
        break;
    }
    if (this->file_ptr==null){
        assert_error(false, "failed to open file");
    }
    this->is_open = true;
}

void file_t::close(){
    if (this->is_open){
        free(this->filename);
        fclose(this->file_ptr);
        this->is_open = false;
    }
}

void file_t::write(const char *format, ...){
    assert_error(this->is_open, "invalid file mode");
    va_list args;
    va_start(args, format);
    vfprintf(this->file_ptr, format, args);
    va_end(args);
}

int_t file_t::read(const char *format, ...){
    assert_error(this->is_open, "can't read from unopened file");
    va_list args;
    va_start(args, format);
    int_t n=vfscanf(this->file_ptr, format, args);
    assert(n>=0);
    va_end(args);
    return feof(this->file_ptr);
}

//

void binary_file_t::open(const char *filename, const char mode){
    assert_error(this->is_open==false, "file is already open");
    this->filename = (char*)calloc(this->max_length, sizeof(char));
    strcpy(this->filename, filename);
    assert_error(mode=='r'||mode=='w'||mode=='a', "invalid file mode");
    this->mode = mode;
    switch (this->mode){
    case 'w':
        this->file_ptr = fopen(this->filename, "wb");
        break;
    case 'r':
        this->file_ptr = fopen(this->filename, "rb");
        break;
    case 'a':
        this->file_ptr = fopen(this->filename, "ab");
        break;
    }
    if (this->file_ptr==null){
        assert_error(false, "failed to open file");
    }
    this->is_open = true;
}
