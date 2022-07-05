


static inline void mul_fp2(quad_t ret, const quad_t a, const quad_t b)
{
 cy_mult_mont(fp_ctx_t *ctx, fp_t *a, fp_t *b);
}

void mul_mont_384x(vec384x ret, const vec384x a, const vec384x b,
                   const vec384 p, limb_t n0)
{
  vec384 aa, bb, cc;

  add_mod_n(aa, a[0], a[1], p, NLIMBS(384));
  add_mod_n(bb, b[0], b[1], p, NLIMBS(384));
  mul_mont_n(bb, bb, aa, p, n0, NLIMBS(384));
  mul_mont_n(aa, a[0], b[0], p, n0, NLIMBS(384));
  mul_mont_n(cc, a[1], b[1], p, n0, NLIMBS(384));
  sub_mod_n(ret[0], aa, cc, p, NLIMBS(384));
  sub_mod_n(ret[1], bb, aa, p, NLIMBS(384));
  sub_mod_n(ret[1], ret[1], cc, p, NLIMBS(384));
}
