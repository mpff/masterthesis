library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
#for basis representation of full proc mean
library(mboost, quietly = TRUE, warn.conflicts = FALSE)
library(shapeboost, quietly = TRUE, warn.conflicts = FALSE)

# load digit 3 example data from shapes package
d3 <- shapes::digit3.dat

## apply relative "arc-length" parametrization
# first move to complex representation for convenience
d3 <- complex(re = d3[,1,], im = d3[,2,]) %>% 
  matrix(nrow =dim(d3)[1])
d3 <-apply(d3, 2, function(x) data.frame(y = x))
d3 <-bind_rows(d3, .id = "id")
# now get parametrization as y(arg)
d3 <- d3 %>% group_by(id )%>% 
  mutate(arg =cumsum(c(0,Mod(diff(y)))) )%>% 
  mutate(arg = arg/max(arg))%>% ungroup()

d3 %>% filter(id%in%1:16) %>% 
  ggplot(aes(Re(y),Im(y)))+ geom_path()+ facet_wrap(~id)

# now represent d3 in a'very long'with y real and
# both dimensions in one column
d3 <-lapply(list(x = Re, y = Im), function(f) mutate(d3, y =f(y))) %>% 
  bind_rows(.id = "dim")

# calculate functional full procrustres mean shape
mean_d3 <- planarshape_full_proc(formula = y^dim~ bbs(arg) | id,
                                 data = d3, 
                                 smoothed.cov = TRUE,
                                 arg.grid.len = 50)

# depict result
mean_d3 %>% as_shape_frame_default() %>% arrange(arg) %>%
  ggplot(aes(real, imaginary))+ geom_path()