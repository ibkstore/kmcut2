fit_ideal_linear_plot<-function(filp.a1,filp.b1,filp.a2,filp.b2,filp.mid_i,filp.max_i)
{
  filp.v = numeric(length = filp.max_i)

  for(filp.i in 1:filp.max_i)
  {
    if(filp.i <= filp.mid_i)
    {
      filp.v[filp.i] = filp.a1 + filp.b1 * filp.i
    }else
    {
      filp.v[filp.i] = filp.a2 + filp.b2 * filp.i
    }
  }
  return(filp.v)
}
