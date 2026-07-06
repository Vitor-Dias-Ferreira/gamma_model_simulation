# ==========================================
# Atividade 3 - Simulação Monte Carlo (Gama)
# ==========================================

# 0. Setup inicial
set.seed(42) # Semente pra bater os mesmos resultados sempre
alpha_0 <- 2 # alpha verdadeiro da pop
beta_0 <- 3  # beta verdadeiro
n_tamanhos <- c(20, 50, 100, 200, 500) # os n's que o pdf pede
R <- 10000   # 10k replicações pra cada n

# 1. Estimador de Momentos (EM) - Forma fechada, mais de boa
calcula_momentos <- function(x) {
  n <- length(x)
  x_barra <- mean(x)
  
  # O R calcula a variância com n-1. Arrumando pro divisor n do trabalho:
  s2 <- var(x) * ((n - 1) / n) 
  
  # Sistema de equações dos momentos resolvido
  alpha_mm <- (x_barra^2) / s2
  beta_mm <- s2 / x_barra
  
  return(c(alpha = alpha_mm, beta = beta_mm))
}

# Testezinho rápido pra ver se não quebra
dados_teste <- rgamma(n = 100, shape = alpha_0, scale = beta_0)
print("Estimativas do Método dos Momentos:")
print(calcula_momentos(dados_teste))


# ==========================================
# 2. Máxima Verossimilhança (Newton-Raphson na unha)
# ==========================================

calcula_emv <- function(x, tol = 1e-8, max_iter = 100) {
  n <- length(x)
  x_barra <- mean(x)
  
  # Essa parte da eq 3 não muda no loop
  constante <- log(x_barra) - mean(log(x))
  
  # Usando os momentos como chute inicial pro algoritmo não se perder
  s2 <- var(x) * ((n - 1) / n)
  alpha_init <- (x_barra^2) / s2
  
  # Passando pro log pra evitar alpha negativo iterando (eta = ln(alpha))
  eta <- log(alpha_init)
  
  for (i in 1:max_iter) {
    alpha_atual <- exp(eta)
    
    # Eq que quero zerar (valeu R por já ter a digamma pronta)
    g_val <- eta - digamma(alpha_atual) - constante
    
    # Derivada pra usar no passo de Newton
    g_deriv <- 1 - trigamma(alpha_atual) * alpha_atual
    
    # Atualiza o eta
    eta_novo <- eta - (g_val / g_deriv)
    
    # Vê se já convergiu
    if (abs(eta_novo - eta) < tol) {
      eta <- eta_novo
      break
    }
    
    eta <- eta_novo
  }
  
  # Voltando pra escala normal antes de retornar
  alpha_mle <- exp(eta)
  beta_mle <- x_barra / alpha_mle
  
  return(c(alpha = alpha_mle, beta = beta_mle))
}

# Testando o EMV pra ver se dá parecido com os Momentos:
print("Estimativas de Máxima Verossimilhança:")
print(calcula_emv(dados_teste))


# ==========================================
# 3. Erros Padrão (Delta, Fisher e Observada)
# ==========================================

calcula_erros_padrao <- function(x, est_mm, est_mle) {
  n <- length(x)
  
  # --- 1. Momentos (MM) ---
  a_mm <- est_mm["alpha"]
  b_mm <- est_mm["beta"]
  
  # Matriz do método delta que tava no pdf
  Sigma_MM <- matrix(c(
    a_mm * (2 + a_mm),        b_mm * (2 + a_mm),
    b_mm * (2 + a_mm),        (b_mm^2 * (2 + a_mm)) / a_mm
  ), nrow = 2, byrow = TRUE)
  
  # Tira a raiz da diagonal pra pegar o erro padrão
  se_mm <- sqrt(diag(Sigma_MM) / n)
  
  # --- 2. EMV Fisher (MLE-IF) ---
  a_mle <- est_mle["alpha"]
  b_mle <- est_mle["beta"]
  
  # Informação de Fisher (Esperada)
  I_F <- matrix(c(
    trigamma(a_mle),   1 / b_mle,
    1 / b_mle,         a_mle / (b_mle^2)
  ), nrow = 2, byrow = TRUE)
  
  # Inverte a matriz multiplicada por n
  se_mle_if <- sqrt(diag(solve(n * I_F)))
  
  # --- 3. EMV Observada (MLE-IO) ---
  # Derivadas segundas pra matriz Hessiana
  H11 <- -n * trigamma(a_mle)
  H12 <- -n / b_mle
  H22 <- (n * a_mle / (b_mle^2)) - (2 * sum(x) / (b_mle^3))
  Hessiana <- matrix(c(H11, H12, H12, H22), nrow = 2, byrow = TRUE)
  
  # Inverte a informação observada
  se_mle_io <- sqrt(diag(solve(-Hessiana)))
  
  return(list(
    SE_MM = setNames(se_mm, c("alpha", "beta")),
    SE_MLE_IF = setNames(se_mle_if, c("alpha", "beta")),
    SE_MLE_IO = setNames(se_mle_io, c("alpha", "beta"))
  ))
}

# Testando se os erros saem bonitinhos:
print("Erros Padrão Calculados:")
est_mm_teste <- calcula_momentos(dados_teste)
est_mle_teste <- calcula_emv(dados_teste)
print(calcula_erros_padrao(dados_teste, est_mm_teste, est_mle_teste))


# ==========================================
# 4. O Loop Pesado (Monte Carlo)
# ==========================================

z_val <- qnorm(0.975) # Z pra 95% de confiança
lista_resultados <- list() # Vai guardar os dfs de cada iteração

# Função pra calcular de uma vez as métricas da rodada pro relatório
calcula_metricas <- function(n_atual, est, se, theta_0, nome_metodo, param) {
  LI <- est - z_val * se
  LS <- est + z_val * se
  
  data.frame(
    n = n_atual,
    Parametro = param,
    Metodo = nome_metodo,
    Media = mean(est),
    Mean_SE = mean(se),
    SD = sd(est),
    Bias = mean(est - theta_0),
    MSE = mean((est - theta_0)^2),
    Largura = mean(LS - LI),
    PC = mean(theta_0 >= LI & theta_0 <= LS), # Proporção que cobriu o valor real
    Alpha1 = mean(theta_0 < LI),              # Ficou fora pela esquerda
    Alpha2 = mean(theta_0 > LS)               # Ficou fora pela direita
  )
}

# Bora rodar tudo
for (n in n_tamanhos) {
  cat(sprintf("Rodando %d simulações para n = %d... (Aguarde)\n", R, n))
  
  # Matrizes R x 2 pra guardar os resultados dessa rodada
  est_mm_mat <- matrix(NA, nrow = R, ncol = 2)
  est_mle_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mm_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mle_if_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mle_io_mat <- matrix(NA, nrow = R, ncol = 2)
  
  for (j in 1:R) {
    # 1. Gera amostra
    x <- rgamma(n, shape = alpha_0, scale = beta_0)
    
    # 2. Salva as ests
    est_mm <- calcula_momentos(x)
    est_mle <- calcula_emv(x)
    
    # 3. Pega os erros
    erros <- calcula_erros_padrao(x, est_mm, est_mle)
    
    # 4. Joga tudo na linha j das matrizes
    est_mm_mat[j, ] <- est_mm
    est_mle_mat[j, ] <- est_mle
    se_mm_mat[j, ] <- erros$SE_MM
    se_mle_if_mat[j, ] <- erros$SE_MLE_IF
    se_mle_io_mat[j, ] <- erros$SE_MLE_IO
  }
  
  # Compilando tudo numa tabelinha pro n atual
  tabela_n <- rbind(
    # --- Alpha ---
    calcula_metricas(n, est_mm_mat[,1], se_mm_mat[,1], alpha_0, "MM", "Alpha"),
    calcula_metricas(n, est_mle_mat[,1], se_mle_if_mat[,1], alpha_0, "MLE-IF", "Alpha"),
    calcula_metricas(n, est_mle_mat[,1], se_mle_io_mat[,1], alpha_0, "MLE-IO", "Alpha"),
    
    # --- Beta ---
    calcula_metricas(n, est_mm_mat[,2], se_mm_mat[,2], beta_0, "MM", "Beta"),
    calcula_metricas(n, est_mle_mat[,2], se_mle_if_mat[,2], beta_0, "MLE-IF", "Beta"),
    calcula_metricas(n, est_mle_mat[,2], se_mle_io_mat[,2], beta_0, "MLE-IO", "Beta")
  )
  
  lista_resultados[[as.character(n)]] <- tabela_n
}

# Empilha as tabelas de cada n num tabelão só
tabela_final <- do.call(rbind, lista_resultados)
rownames(tabela_final) <- NULL

print("Simulação Concluída! Aqui estão as primeiras linhas da tabela consolidada:")
print(head(tabela_final))


# ==========================================
# 5. Exportação pro LaTeX e Gráficos
# ==========================================
library(ggplot2)

# Pasta pra não virar zona no diretório principal
pasta_resultados <- "resultados_gama"
if(!dir.exists(pasta_resultados)) dir.create(pasta_resultados)

# O CSV que vai salvar a vida na hora de montar a tabela no LaTeX
caminho_tabela <- file.path(pasta_resultados, "tabela_completa.csv")
write.csv(tabela_final, caminho_tabela, row.names = FALSE)
cat("Tabela salva em:", caminho_tabela, "\n")

# Pegando só MM e MLE-IF do Alpha, senão o plot vira um carnaval
dados_grafico <- subset(tabela_final, Parametro == "Alpha" & Metodo %in% c("MM", "MLE-IF"))

# Plot 1: MSE provando a eficiência do MLE
g_mse <- ggplot(dados_grafico, aes(x = n, y = MSE, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Erro Quadrático Médio (MSE) vs n",
       subtitle = "Estimador do parâmetro de forma (Alpha)",
       y = "MSE", x = "Tamanho da Amostra (n)") +
  scale_color_brewer(palette = "Set1")

# Plot 2: A cobertura bate os 95% assintóticos?
g_pc <- ggplot(dados_grafico, aes(x = n, y = PC, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + 
  theme_minimal() +
  labs(title = "Probabilidade de Cobertura vs n",
       subtitle = "Estimador Alpha (Confiança de 95%)",
       y = "Cobertura", x = "Tamanho da Amostra (n)") +
  scale_color_brewer(palette = "Set1")

# Plot 3: IC ficando mais fino com n maior
g_largura <- ggplot(dados_grafico, aes(x = n, y = Largura, color = Metodo, group = Metodo)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Largura Média dos Intervalos vs n",
       subtitle = "Estimador Alpha",
       y = "Largura do IC", x = "Tamanho da Amostra (n)") +
  scale_color_brewer(palette = "Set1")

# Joga tudo pra png direto na pasta pra importar no .tex
ggsave(file.path(pasta_resultados, "grafico_mse.png"), plot = g_mse, width = 8, height = 5)
ggsave(file.path(pasta_resultados, "grafico_pc.png"), plot = g_pc, width = 8, height = 5)
ggsave(file.path(pasta_resultados, "grafico_largura.png"), plot = g_largura, width = 8, height = 5)

cat("Gráficos gerados na pasta:", pasta_resultados, "\n")