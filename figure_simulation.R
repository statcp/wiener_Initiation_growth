library(openxlsx)
library(ggplot2)
library(tidyverse)
sum_dat = list()
for(i in 1:2){
  sum_dat[[i]] = read.xlsx("Data_summary.xlsx",sheet = i)
}

library(ggsci)
library(viridis)
library(latex2exp)
library(grid)
library(gridExtra)
# install.packages("patchwork")
library(patchwork)

f_names <- list('mu0' = TeX(c("$\\mu_{0}$")), 'sig0' = TeX(c("$\\sigma_{0}^2$")),
                'a' = TeX(c("$a$")), 'b' = TeX(c("$b$")),
                'sigb' = TeX(c("$\\sigma_{B}^2$")), 're' = TeX(c("$R(t_{0})=0.8$")))
f_labeller <- function(variable, value){return(f_names[value])}

# =============== Figure of coverage probabilities for parameters ===============

data.frame(sum_dat[[1]]) %>% filter(content == "cp",level == '95') -> t1cp_dat #线性 p=2
colnames(t1cp_dat) = c("n","m","Method","CP","level","mu0","sig0","a","b","sigb","re")
t1cp_dat %>% pivot_longer(mu0:re,values_to = "value", names_to = "name") -> t1cp_dat_new
t1cp_dat_new$name = factor(t1cp_dat_new$name,levels = c("mu0","sig0","a","b","sigb","re"))


fig_cp_pa = t1cp_dat_new %>% 
  ggplot(aes(factor(n), value,color = factor(Method)))  + 
  geom_point(size=2)+ 
  geom_line(aes(group = Method),linetype = 2) +
  facet_wrap(vars(factor(name)),scales = "free", labeller = f_labeller, nrow=2) +
  theme_bw() + 
  labs(color = "Method")+
  scale_color_manual(name = "Method", values = c("blue","red")) +
  theme(legend.position = "right") + xlab("n") + ylab(TeX(r'(CP)'))+ylim(0.93,0.97)

fig_cp_pa


# =============== Figure of average lengths for parameters ===============

data.frame(sum_dat[[1]]) %>% filter(content == "len",level == '95') -> t1len_dat #线性 p=2
colnames(t1len_dat) = c("n","m","Method","Len","level","mu0","sig0","a","b","sigb","re")
t1len_dat %>% pivot_longer(mu0:re,values_to = "value", names_to = "name") -> t1len_dat_new
t1len_dat_new$name = factor(t1len_dat_new$name,levels = c("mu0","sig0","a","b","sigb","re"))
fig_len_pa = t1len_dat_new %>% 
  ggplot(aes(factor(n), value,color = factor(Method)))  + 
  geom_point(size=2)+ 
  geom_line(aes(group = Method),linetype = 2) +
  facet_wrap(vars(factor(name)),scales = "free_y", labeller = f_labeller, nrow=2) +
  theme_bw() + 
  labs(color = "Method")+
  scale_color_manual(name = "Method", values = c("blue","red")) +
  theme(legend.position = "right") + xlab("n") + ylab(TeX(r'(Length)'))
fig_len_pa




# =============== Figure of coverage probabilities and average lengths for RULs  ===============
data.frame(sum_dat[[2]]) %>% filter(content == "cp",level == '95') -> t1rulcp_dat #线性 p=2
colnames(t1rulcp_dat) = c("n","m","Method","CP","level","yk1","yk2","yk3")
f_names <- list('yk1' = TeX(c("$y_k = \\omega/3 $")), 'yk2' = TeX(c("$y_k =  \\omega/2$")),
                'yk3' = TeX(c("$y_k = 2\\omega/3$")))
t1rulcp_dat %>% pivot_longer(yk1:yk3,values_to = "value", names_to = "name") -> t1rulcp_dat_new

data.frame(sum_dat[[2]]) %>% filter(content == "len",level == '95') -> t1rullen_dat #线性 p=2
colnames(t1rullen_dat) = c("n","m","Method","Len","level","yk1","yk2","yk3")
t1rullen_dat %>% pivot_longer(yk1:yk3,values_to = "value", names_to = "name") -> t1rullen_dat_new


fig_cp_rul = t1rulcp_dat_new %>% 
  ggplot(aes(factor(n), value,color = factor(Method)))  + 
  geom_point(size=2)+ 
  geom_line(aes(group = Method),linetype = 2) +
  facet_wrap(vars(factor(name)),scales = "fixed", labeller = f_labeller, nrow=1) +
  theme_bw() + 
  labs(color = "Method")+
  scale_color_manual(name = "Method", values = c("blue","red")) +
  theme(legend.position = "right") +  xlab(NULL) +ylab(TeX(r'(CP)'))+ylim(0.91,0.98)
fig_cp_rul



fig_len_rul = t1rullen_dat_new %>% 
  ggplot(aes(factor(n), value,color = factor(Method)))  + 
  geom_point(size=2)+ 
  geom_line(aes(group = Method),linetype = 2) +
  facet_wrap(vars(factor(name)),scales = "fixed", labeller = f_labeller, nrow=1) +
  theme_bw() + 
  labs(color = "Method")+
  scale_color_manual(name = "Method", values = c("blue","red")) +
  theme(legend.position = "right") + xlab("n") + ylab(TeX(r'(Length)'))

fig_len_rul

p1rul = fig_len_cp /fig_len_rul + plot_layout(ncol = 1, guides = "collect")
p1rul
