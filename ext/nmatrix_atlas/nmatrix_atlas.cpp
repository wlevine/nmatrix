#include <ruby.h>

VALUE cNMatrix;

typedef VALUE (*METHOD)(...);

extern "C" {

static VALUE nm_test_atlas(VALUE self) {
  return INT2FIX(0);
}

void Init_nmatrix_atlas() {
  cNMatrix = rb_define_class("NMatrix", rb_cObject);

	rb_define_method(cNMatrix, "test_atlas", (METHOD)nm_test_atlas, 0);
}

}
