# BART
Bayesian additive tree 코드
######TREE 관련 함수######
#Node 정의
create_node <- function(is_terminal = TRUE,
                        split_variable = NULL,
                        split_value = NULL,
                        left = NULL,
                        right = NULL,
                        data_indices = integer(0),
                        depth = 0,
                        node_id = 1,
                        parent_id = NA) {
  node <- list(is_terminal = is_terminal,
               split_variable = split_variable,
               split_value = split_value,
               left = left,
               right = right,
               data_indices = data_indices,
               depth = depth,
               node_id = as.numeric(node_id),
               parent_id = as.numeric(parent_id)
  )
  class(node) <- "BCARTNode"
  return(node)
}

#Tree 정의
create_tree <- function(root_node) {
  tree <- list(root = root_node,
               nodes = list(root_node)
  )
  class(tree) <- "BCARTree"
  return(tree)
}

#초기 트리 생성 함수
create_initial_tree <- function(X, N, v) {
  root <- create_node(
    is_terminal = TRUE,
    data_indices = 1:nrow(X),
    depth = 0,
    node_id = 1
  )
  
  tree <- create_tree(root)
  return(tree)
}


#GROW 함수 정의
#split variable 선택
get_split_candidates <- function(X, indices) {
  valid_vars <- c()
  for (j in seq_along(X)) {
    x_sub <- X[[j]][indices]
    if (length(unique(x_sub)) > 1) {
      valid_vars <- c(valid_vars, j)
    }
  }
  return(valid_vars)
}

#split value 선택
choose_split_value <- function(x_sub) {
  if (is.factor(x_sub)) {
    lvls <- unique(x_sub)
    k <- sample(1:(length(lvls)-1), 1)
    subset <- sample(lvls, k)
    return(subset)  # return level subset
  } else {
    vals <- sort(unique(x_sub))
    midpoints <- (vals[-length(vals)] + vals[-1]) / 2
    return(sample(midpoints, 1))
  }
}


#GROW 함수
grow_tree <- function(tree, X) {
  terminal_nodes <- Filter(function(n) n$is_terminal, tree$nodes)
  if (length(terminal_nodes) == 0) return(tree)
  
  # 무작위로 하나 선택
  selected_node <- sample(terminal_nodes, 1)[[1]]
  data_indices <- selected_node$data_indices
  
  # 가능한 split 변수 확인
  split_vars <- get_split_candidates(X, data_indices)
  if (length(split_vars) == 0) return(tree)
  
  # split 변수 및 값 선택
  split_var <- sample(split_vars, 1)
  x_sub <- X[[split_var]][data_indices]
  split_val <- choose_split_value(x_sub)
  
  # split 적용
  if (is.factor(X[[split_var]])) {
    left_idx <- data_indices[X[[split_var]][data_indices] %in% split_val]
    right_idx <- data_indices[!X[[split_var]][data_indices] %in% split_val]
  } else {
    left_idx <- data_indices[X[[split_var]][data_indices] <= split_val]
    right_idx <- data_indices[X[[split_var]][data_indices] > split_val]
  }
  
  # 최소 분할 크기 조건 확인
  if (length(left_idx) < 2 || length(right_idx) < 2) {
    return(tree)
  }
  
  # 새로운 노드 생성
  max_id <- max(sapply(tree$nodes, function(n) n$node_id))
  left_node <- create_node(data_indices = left_idx,
                           depth = selected_node$depth + 1,
                           node_id = max_id + 1,
                           parent_id = selected_node$node_id)
  
  right_node <- create_node(data_indices = right_idx,
                            depth = selected_node$depth + 1,
                            node_id = max_id + 2,
                            parent_id = selected_node$node_id)
  
  # 선택된 노드 갱신
  selected_node$is_terminal <- FALSE
  selected_node$split_variable <- split_var
  selected_node$split_value <- split_val
  selected_node$left <- left_node
  selected_node$right <- right_node
  
  # tree 구조에 반영
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == selected_node$node_id))]] <- selected_node
  tree$nodes <- c(tree$nodes, list(left_node, right_node))
  
  if (selected_node$node_id == tree$root$node_id) {
    tree$root <- selected_node
  }
  
  # 연결 관계 재설정
  tree <- update_tree_structure(tree)
  
  # Split 내용 출력
  cond_str <- if (is.numeric(split_val)) {
    paste0("X", split_var, " <= ", split_val)
  } else {
    paste0("X", split_var, " in {", paste(split_val, collapse = ","), "}")
  }
  return(tree)
}

#PRUNE 함수
#Terminal node인지 확인
is_prunable <- function(node) {
  !node$is_terminal &&
    !is.null(node$left) && !is.null(node$right) &&
    node$left$is_terminal && node$right$is_terminal
}

#PRUNE
prune_tree <- function(tree) {
  # prune 가능한 부모 노드 리스트
  prunable_nodes <- Filter(is_prunable, tree$nodes)
  if (length(prunable_nodes) == 0) return(tree)
  
  # 무작위로 하나 선택
  parent_node <- sample(prunable_nodes, 1)[[1]]
  left_id <- parent_node$left$node_id
  right_id <- parent_node$right$node_id
  
  # 해당 자식 노드 제거
  tree$nodes <- Filter(function(n) !(n$node_id %in% c(left_id, right_id)), tree$nodes)
  
  # 부모 노드를 terminal로 전환
  parent_node$is_terminal <- TRUE
  parent_node$split_variable <- NULL
  parent_node$split_value <- NULL
  parent_node$left <- NULL
  parent_node$right <- NULL
  
  # 노드 정보 갱신
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == parent_node$node_id))]] <- parent_node
  
  tree <- update_tree_structure(tree)
  
  return(tree)
}


#CHANGE 함수
#내부 노드 중 하나 선택
get_internal_nodes <- function(tree) {
  Filter(function(n) !n$is_terminal, tree$nodes)
}

#CHANGE 1 : 분할 값만 변경
change1_tree <- function(tree, X) {
  internal_nodes <- get_internal_nodes(tree)
  if (length(internal_nodes) == 0) return(tree)
  
  node <- sample(internal_nodes, 1)[[1]]
  var_idx <- node$split_variable
  if (is.null(var_idx)) return(tree)
  
  data_idx <- node$data_indices
  x_vals <- X[data_idx, var_idx]
  new_split <- choose_split_value(x_vals)
  
  # 동일하면 변경 의미 없음
  if (is.factor(x_vals)) {
    if (setequal(new_split, node$split_value)) return(tree)
  } else {
    if (new_split == node$split_value) return(tree)
  }
  
  node$split_value <- new_split
  
  # 기존 자식 노드와 구조는 유지, 단 split value만 변경
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == node$node_id))]] <- node
  
  tree <- update_tree_structure(tree)
  
  return(tree)
}

#CHANGE 2 : 변수와 분할 값 모두 변경
change2_tree <- function(tree, X) {
  internal_nodes <- get_internal_nodes(tree)
  if (length(internal_nodes) == 0) return(tree)
  
  node <- sample(internal_nodes, 1)[[1]]
  data_idx <- node$data_indices
  
  # 가능한 변수 찾기
  split_vars <- get_split_candidates(X, data_idx)
  if (length(split_vars) == 0) return(tree)
  
  var_idx <- sample(split_vars, 1)
  new_split <- choose_split_value(X[[var_idx]][data_idx])
  
  same_split <- FALSE
  if (is.factor(X[[var_idx]])) {
    same_split <- setequal(new_split, node$split_value)
  } else {
    same_split <- new_split == node$split_value
  }
  
  if (var_idx == node$split_variable && same_split) return(tree)
  
  
  node$split_variable <- var_idx
  node$split_value <- new_split
  
  # 기존 자식 구조 유지, 단 split rule만 변경
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == node$node_id))]] <- node
  
  tree <- update_tree_structure(tree)
  
  return(tree)
}

# SWAP 
swap_tree <- function(tree) {
  internal_nodes <- get_internal_nodes(tree)
  
  swap_candidates <- list()
  
  for (parent in internal_nodes) {
    if (!is.null(parent$left) && !parent$left$is_terminal) {
      # 분할 기준이 다른 경우만 후보로 추가
      if (!(parent$split_variable == parent$left$split_variable &&
            parent$split_value == parent$left$split_value)) {
        swap_candidates <- append(swap_candidates, list(list(parent, parent$left)))
      }
    }
    if (!is.null(parent$right) && !parent$right$is_terminal) {
      if (!(parent$split_variable == parent$right$split_variable &&
            parent$split_value == parent$right$split_value)) {
        swap_candidates <- append(swap_candidates, list(list(parent, parent$right)))
      }
    }
  }
  
  if (length(swap_candidates) == 0) {
    return(tree)
  }
  
  # 유효한 쌍 중에서 무작위 선택
  chosen_pair <- sample(swap_candidates, 1)[[1]]
  parent_node <- chosen_pair[[1]]
  child_node <- chosen_pair[[2]]
  
  # 디버깅 로그
  cat("SWAP executed:\n")
  cat(" Before → Parent node", parent_node$node_id, ": X", parent_node$split_variable, "≤", parent_node$split_value, "\n")
  cat("           Child node", child_node$node_id, ": X", child_node$split_variable, "≤", child_node$split_value, "\n")
  
  # 교환
  tmp_var <- parent_node$split_variable
  tmp_val <- parent_node$split_value
  parent_node$split_variable <- child_node$split_variable
  parent_node$split_value <- child_node$split_value
  child_node$split_variable <- tmp_var
  child_node$split_value <- tmp_val
  
  cat(" After  → Parent node", parent_node$node_id, ": X", parent_node$split_variable, "≤", parent_node$split_value, "\n")
  cat("           Child node", child_node$node_id, ": X", child_node$split_variable, "≤", child_node$split_value, "\n")
  
  # 갱신
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == parent_node$node_id))]] <- parent_node
  tree$nodes[[which(sapply(tree$nodes, function(n) n$node_id == child_node$node_id))]] <- child_node
  
  tree <- update_tree_structure(tree)
  return(tree)
}


#트리 업데이트
update_tree_structure <- function(tree) {
  id_map <- setNames(tree$nodes, as.character(sapply(tree$nodes, function(n) n$node_id)))
  
  reconnect <- function(node) {
    if (!node$is_terminal) {
      left_id <- as.character(node$left$node_id)
      right_id <- as.character(node$right$node_id)
      
      if (is.null(id_map[[left_id]])) {
        cat("⚠️ Left child ID", left_id, "not found in id_map\n")
        return(node)
      }
      if (is.null(id_map[[right_id]])) {
        cat("⚠️ Right child ID", right_id, "not found in id_map\n")
        return(node)
      }
      
      node$left <- reconnect(id_map[[left_id]])
      node$right <- reconnect(id_map[[right_id]])
    }
    return(node)
  }
  
  tree$root <- reconnect(id_map[["1"]])
  return(tree)
}

print_all_nodes <- function(tree) {
  for (node in tree$nodes) {
    cat("Node ID:", node$node_id,
        if (node$is_terminal) "[Terminal]" else "[Internal]",
        "| Depth:", node$depth,
        "| Data size:", length(node$data_indices), "\n")
    
    if (!node$is_terminal) {
      cat("  Split: X", node$split_variable, " <= ", node$split_value, "\n")
      cat("  Left child ID:", node$left$node_id, "\n")
      cat("  Right child ID:", node$right$node_id, "\n")
    }
  }
}



#트리 구조 출력 함수
print_tree <- function(node, indent = "") {
  cat(indent, "Node ID:", node$node_id,
      if (node$is_terminal) "[Terminal]" else "[Internal]",
      "| Data size:", length(node$data_indices), "\n")
  if (!node$is_terminal) {
    cat(indent, " Split: X", node$split_variable, " <= ", node$split_value, "\n")
    print_tree(node$left, paste0(indent, "  "))
    print_tree(node$right, paste0(indent, "  "))
  }
}

######ZIP-BCART 관련 함수######
#log-prior
compute_tree_prior <- function(tree, alpha = 0.95, beta = 2.0) {
  log_prior <- 0  # 로그 확률로 계산 (언더플로 방지)
  
  for (node in tree$nodes) {
    depth <- node$depth
    if (node$is_terminal) {
      # 단말 노드일 경우: 분할 안 될 확률
      log_prior <- log_prior + log(1 - alpha * (1 + depth)^(-beta))
    } else {
      # 내부 노드일 경우: 분할될 확률
      log_prior <- log_prior + log(alpha * (1 + depth)^(-beta))
    }
  }
  
  return(log_prior)  # 로그 prior 값을 반환
}

#log-integrated augmented likelihood for terminal node t
compute_integrated_aug_likelihood_zip1 <- function(node, N, v, delta, phi,
                                                   alpha1 = 1, beta1 = 1,
                                                   alpha2 = 1, beta2 = 1) {
  idx <- node$data_indices
  if (length(idx) == 0) return(-Inf)  # 공백 노드는 로그 우도 -Inf 처리
  
  N_t <- N[idx]
  v_t <- v[idx]
  delta_t <- delta[idx]
  phi_t <- phi[idx]
  
  # 항별 계산
  sum_delta <- sum(delta_t)
  sum_phi <- sum(phi_t)
  sum_dN <- sum(delta_t * N_t)
  sum_dv <- sum(delta_t * v_t)
  
  # log-우도 계산: 상수 × 지수 × 감마 분포 적분 결과
  log_const <- alpha1 * log(beta1) - lgamma(alpha1) +
    alpha2 * log(beta2) - lgamma(alpha2)
  
  log_prod_part <- sum(-phi_t) +
    sum(delta_t * N_t * log(v_t + 1e-10)) -
    sum(delta_t * lgamma(N_t + 1))
  
  log_gamma1 <- lgamma(sum_delta + alpha1) -
    (sum_delta + alpha1) * log(sum_phi + beta1)
  
  log_gamma2 <- lgamma(sum_dN + alpha2) -
    (sum_dN + alpha2) * log(sum_dv + beta2)
  
  # 최종 log likelihood
  log_likelihood <- log_const + log_prod_part + log_gamma1 + log_gamma2
  
  return(log_likelihood)
}

#Integrated augmented likelihood for the tree
compute_tree_likelihood_zip1 <- function(tree, N, v, delta, phi,
                                         alpha1 = 1, beta1 = 1,
                                         alpha2 = 1, beta2 = 1) {
  log_likelihood <- 0
  
  for (node in tree$nodes) {
    if (node$is_terminal) {
      log_likelihood <- log_likelihood + compute_integrated_aug_likelihood_zip1(
        node = node,
        N = N,
        v = v,
        delta = delta,
        phi = phi,
        alpha1 = alpha1,
        beta1 = beta1,
        alpha2 = alpha2,
        beta2 = beta2
      )
    }
  }
  
  return(log_likelihood)  # 로그 우도 반환
}

#Posterior
compute_tree_posterior_zip1 <- function(tree, N, v, delta, phi,
                                        alpha1 = 1, beta1 = 1,
                                        alpha2 = 1, beta2 = 1,
                                        alpha_tree = 0.95, beta_tree = 2.0) {
  log_prior <- compute_tree_prior(tree, alpha = alpha_tree, beta = beta_tree)
  log_likelihood <- compute_tree_likelihood_zip1(
    tree = tree,
    N = N,
    v = v,
    delta = delta,
    phi = phi,
    alpha1 = alpha1,
    beta1 = beta1,
    alpha2 = alpha2,
    beta2 = beta2
  )
  
  log_posterior <- log_prior + log_likelihood
  return(log_posterior)
}

#MH 알고리즘
mcmc_zip1_tree_with_latent_and_parameters <- function(X, N, v,
                                                      n_iter = 100,
                                                      alpha1 = 1, beta1 = 1,
                                                      alpha2 = 1, beta2 = 1,
                                                      alpha_tree = 0.95, beta_tree = 2.0) {
  n <- length(N)
  delta <- rep(1, n)
  phi <- rexp(n, rate = 1 + 1)
  
  root_node <- create_node(data_indices = 1:n)
  tree <- create_tree(root_node)
  
  posterior_trace <- numeric(n_iter)
  delta_chain <- matrix(NA, nrow = n_iter, ncol = n)
  phi_chain <- matrix(NA, nrow = n_iter, ncol = n)
  
  mu_chain <- vector("list", n_iter)
  lambda_chain <- vector("list", n_iter)
  node_id_chain <- vector("list", n_iter)
  
  best_tree <- tree
  best_posterior <- -Inf
  
  moves <- c("GROW", "PRUNE", "CHANGE1", "CHANGE2","SWAP")
  
  for (iter in 1:n_iter) {
    latent <- update_latent_variables(tree, N, v,
                                      alpha1 = alpha1, beta1 = beta1,
                                      alpha2 = alpha2, beta2 = beta2)
    delta <- latent$delta
    phi <- latent$phi
    
    delta_chain[iter, ] <- delta
    phi_chain[iter, ] <- phi
    
    param <- sample_node_parameters(tree, N, v, delta, phi,
                                    alpha1, beta1, alpha2, beta2)
    
    # terminal node에 mu/lambda 저장
    for (k in seq_along(param$node_id)) {
      idx <- which(sapply(tree$nodes, function(n) n$node_id) == param$node_id[k])
      if (length(idx) == 1) {
        tree$nodes[[idx]]$mu <- param$mu[k]
        tree$nodes[[idx]]$lambda <- param$lambda[k]
      }
    }
    
    mu_chain[[iter]] <- param$mu
    lambda_chain[[iter]] <- param$lambda
    node_id_chain[[iter]] <- param$node_id
    
    current_post <- compute_tree_posterior_zip1(tree, N, v, delta, phi,
                                                alpha1, beta1, alpha2, beta2,
                                                alpha_tree, beta_tree)
    
    move <- sample(moves, 1,prob = c(0.3,0.3,0.15,0.15,0.1))
    proposed_tree <- switch(move,
                            GROW = grow_tree(tree, X),
                            PRUNE = prune_tree(tree),
                            CHANGE1 = change1_tree(tree, X),
                            CHANGE2 = change2_tree(tree, X),
                            SWAP = swap_tree(tree))
    
    proposed_post <- compute_tree_posterior_zip1(proposed_tree, N, v, delta, phi,
                                                 alpha1, beta1, alpha2, beta2,
                                                 alpha_tree, beta_tree)
    
    log_accept_ratio <- proposed_post - current_post
    if (log(runif(1)) < log_accept_ratio) {
      cat("accepted\n")
      tree <- proposed_tree
      current_post <- proposed_post
    } else { cat("rejected\n")
    }
    
    if (current_post > best_posterior) {
      best_tree <- tree
      best_posterior <- current_post
    }
    
    posterior_trace[iter] <- current_post
  }
  
  #best_tree에 parameter 재할당
  final_param <- sample_node_parameters(best_tree, N, v, delta, phi,
                                        alpha1, beta1, alpha2, beta2)
  
  for (k in seq_along(final_param$node_id)) {
    idx <- which(sapply(best_tree$nodes, function(n) n$node_id) == final_param$node_id[k])
    if (length(idx) == 1) {
      best_tree$nodes[[idx]]$mu <- final_param$mu[k]
      best_tree$nodes[[idx]]$lambda <- final_param$lambda[k]
    }
  }
  
  return(list(
    tree = best_tree,
    trace = posterior_trace,
    delta_chain = delta_chain,
    phi_chain = phi_chain,
    mu_chain = mu_chain,
    lambda_chain = lambda_chain,
    node_id_chain = node_id_chain
  ))
}



#잠재변수 업데이트
update_latent_variables <- function(tree, N, v, alpha1 = 1, beta1 = 1, alpha2 = 1, beta2 = 1) {
  n <- length(N)
  delta_new <- numeric(n)
  phi_new <- numeric(n)
  
  for (node in tree$nodes) {
    if (!node$is_terminal) next
    idx <- node$data_indices
    N_t <- N[idx]
    v_t <- v[idx]
    
    # sufficient statistics
    # 초기값 (전 iteration 기준 delta=1이라고 가정)
    sum_phi <- length(idx) / (1 + 1)  # 1 / E[1 + mu]
    sum_dN <- sum(N_t)
    sum_dv <- sum(v_t)
    sum_d <- length(idx)
    
    # MAP 근사
    mu_hat <- (sum_d + alpha1) / (sum_phi + beta1)
    lambda_hat <- (sum_dN + alpha2) / (sum_dv + beta2)
    
    for (i in seq_along(idx)) {
      ii <- idx[i]
      
      # phi_i ~ Exp(1 + mu_hat)
      phi_new[ii] <- rexp(1, rate = 1 + mu_hat)
      
      # delta_i | N_i
      if (N[ii] == 0) {
        p1 <- mu_hat * exp(-lambda_hat * v[ii])
        p0 <- 1
        prob_1 <- p1 / (p1 + p0 + 1e-10)
        delta_new[ii] <- rbinom(1, 1, prob_1)
      } else {
        delta_new[ii] <- 1  # zero-inflated part는 0만 가능
      }
    }
  }
  
  return(list(delta = delta_new, phi = phi_new))
}

#mu, lambda sampling
sample_node_parameters <- function(tree, N, v, delta, phi,
                                   alpha1, beta1, alpha2, beta2) {
  mu <- rep(NA, length(N))
  lambda <- rep(NA, length(N))
  node_id <- rep(NA, length(N))
  
  for (i in seq_along(tree$nodes)) {
    node <- tree$nodes[[i]]
    
    if (node$is_terminal) {
      indices <- node$data_indices
      
      # 잠재 변수 부분합 계산
      sum_delta <- sum(delta[indices])
      sum_phi   <- sum(phi[indices])
      sum_N     <- sum(N[indices])
      sum_v     <- sum(v[indices])
      
      # mu 샘플링
      mu_sample <- rgamma(1, shape = alpha1 + sum_delta,
                          rate  = beta1 + sum_phi)
      
      # lambda 샘플링
      lambda_sample <- rgamma(1, shape = alpha2 + sum_N,
                              rate  = beta2 + sum_v)
      
      # ✅ 노드에 직접 저장
      node$mu <- mu_sample
      node$lambda <- lambda_sample
      
      # ✅ 트리 객체에도 반영
      tree$nodes[[i]] <- node
      
      # 관측치별 벡터에도 복사
      mu[indices]     <- mu_sample
      lambda[indices] <- lambda_sample
      node_id[indices] <- node$node_id
    }
  }
  
  return(list(mu = mu, lambda = lambda, node_id = node_id, tree = tree))
}


#best_tree의 사후 파라미터
print_terminal_node_parameters <- function(tree) {
  for (node in tree$nodes) {
    if (!is.null(node$mu) && node$is_terminal) {
      cat(sprintf("Node ID %d | Size = %d | mu = %.3f | lambda = %.3f | mu*lambda = %.3f\n",
                  node$node_id,
                  length(node$data_indices),
                  node$mu,
                  node$lambda,
                  node$mu * node$lambda))
    }
  }
}

#tree에서 파라미터(mu,lambda) 가져오기
compute_mu_from_tree <- function(tree, X) {
  mu <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    node <- find_terminal_node(tree, X[i, ,drop = FALSE])
    if (is.null(node$mu)) stop(paste("mu not found for observation", i))
    mu[i] <- node$mu
  }
  return(mu)
}

compute_lambda_from_tree <- function(tree, X) {
  n <- nrow(X)
  lambda <- numeric(n)
  
  for (i in 1:n) {
    node <- find_terminal_node(tree, X[i, , drop = FALSE])  # 여기서 통일
    if (is.null(node$lambda)) stop(paste("lambda not found for observation", i))
    lambda[i] <- node$lambda
  }
  return(lambda)
}

#관측치가 속한 terminal node 찾기
find_terminal_node <- function(tree, x) {
  node <- tree$root
  
  while (!node$is_terminal) {
    var <- node$split_variable
    val <- node$split_value
    xval <- x[[var]]
    
    go_left <- if (is.factor(xval)) {
      xval %in% val
    } else {
      xval <= val
    }
    
    child <- if (go_left) node$left else node$right
    child_id <- if (inherits(child, "BCARTNode")) child$node_id else child
    
    # 🔧 여기 수정
    node_ids <- sapply(tree$nodes, function(n) as.numeric(n$node_id))
    idx <- which(node_ids == as.numeric(child_id))  # ✅ 안전하게 변환
    
    if (length(idx) == 0) stop(paste("Node ID", child_id, "not found"))
    node <- tree$nodes[[idx]]
  }
  
  return(node)
}

#zip_bcart 예측 함수
predict_zipbcart <- function(tree, X_test, v_test) {
  n <- nrow(X_test)
  mu <- numeric(n)
  lambda <- numeric(n)
  
  for (i in 1:n) {
    node <- find_terminal_node(tree, X_test[i, , drop = FALSE])
    
    if (is.null(node$mu) || is.null(node$lambda)) {
      stop(paste("Missing mu/lambda in terminal node for obs", i))
    }
    
    mu[i] <- node$mu
    lambda[i] <- node$lambda
  }
  
  y_pred <- (mu/(1+mu)) * lambda * v_test
  return(y_pred)
}

#####ZIP-BART 관련 함수######
# ZIP-BART prior log-density for a single tree
compute_zip_bart_tree_prior <- function(tree, alpha1, beta1, alpha2, beta2, alpha_tree = 0.95, beta_tree = 1.0) {
  log_prior <- 0
  for (node in tree$nodes) {
    depth <- node$depth
    if (node$is_terminal) {
      log_prior <- log_prior + log(1 - alpha_tree * (1 + depth)^(-beta_tree))
      if (!is.null(node$mu) && !is.null(node$lambda)) {
        log_prior <- log_prior + dgamma(node$mu, shape = alpha1, rate = beta1, log = TRUE) +
          dgamma(node$lambda, shape = alpha2, rate = beta2, log = TRUE)
      }
    } else {
      log_prior <- log_prior + log(alpha_tree * (1 + depth)^(-beta_tree))
    }
  }
  return(log_prior)
}

#residual 기반 likelihood
compute_tree_likelihood_zip_residual <- function(tree,
                                                 mu_residual,
                                                 lambda_residual,
                                                 v,
                                                 alpha1, beta1,
                                                 alpha2, beta2) {
  log_likelihood <- 0
  
  for (node in tree$nodes) {
    if (node$is_terminal) {
      indices <- node$data_indices
      
      # 잔차의 합
      sum_mu <- sum(mu_residual[indices])
      sum_lambda <- sum(lambda_residual[indices])
      sum_v <- sum(v[indices])
      n <- length(indices)
      
      # mu: Gamma 마르지널 로그 가능도 (conjugate case)
      ll_mu <- lgamma(alpha1 + sum_mu) - (alpha1 + sum_mu) * log(beta1 + n) +
        alpha1 * log(beta1) - lgamma(alpha1)
      
      # lambda: Gamma 마르지널 로그 가능도
      ll_lambda <- lgamma(alpha2 + sum_lambda) - (alpha2 + sum_lambda) * log(beta2 + sum_v) +
        alpha2 * log(beta2) - lgamma(alpha2)
      
      log_likelihood <- log_likelihood + ll_mu + ll_lambda
    }
  }
  
  return(log_likelihood)
}

# ZIP-BART log posterior 계산 (prior + likelihood)
compute_tree_posterior_zipbart <- function(tree,
                                           N, v, delta, phi,
                                           mu_residual, lambda_residual,
                                           alpha1 = 1, beta1 = 5, alpha2 = 1, beta2 = 1,
                                           alpha_tree = 0.95, beta_tree = 1.0) {
  # Tree 구조 prior (그대로 유지)
  log_prior <- compute_zip_bart_tree_prior(tree, alpha1, beta1, alpha2, beta2, alpha_tree, beta_tree)
  
  # Residual 기반 log-likelihood 계산
  log_likelihood <- compute_tree_likelihood_zip_residual(
    tree,
    mu_residual = mu_residual,
    lambda_residual = lambda_residual,
    v = v,
    alpha1 = alpha1, beta1 = beta1,
    alpha2 = alpha2, beta2 = beta2
  )
  
  return(log_prior + log_likelihood)
}

# 각 트리에 대해 업데이트 함수 (ZIP likelihood 기반)
update_zip_bart_tree <- function(tree, X, N, v, delta, phi,
                                 mu_minus_j, lambda_minus_j,
                                 mu_residual, lambda_residual,
                                 alpha1, beta1, alpha2, beta2,
                                 alpha_tree, beta_tree) {
  
  # 현재 posterior 계산 (residual 기반으로 계산해야 함)
  current_post <- compute_tree_posterior_zipbart(tree, N, v, delta, phi,
                                                 mu_residual, lambda_residual,
                                                 alpha1, beta1, alpha2, beta2,
                                                 alpha_tree, beta_tree)
  
  move <- sample(c("GROW", "PRUNE", "CHANGE1", "CHANGE2", "SWAP"), 1, prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
  proposed_tree <- switch(move,
                          GROW = grow_tree(tree, X),
                          PRUNE = prune_tree(tree),
                          CHANGE1 = change1_tree(tree, X),
                          CHANGE2 = change2_tree(tree, X),
                          SWAP = swap_tree(tree))
  
  proposed_post <- compute_tree_posterior_zipbart(proposed_tree, N, v, delta, phi,
                                                  mu_residual, lambda_residual,
                                                  alpha1, beta1, alpha2, beta2,
                                                  alpha_tree, beta_tree)
  
  if (log(runif(1)) < (proposed_post - current_post)) {
    tree <- proposed_tree
  }
  
  param <- sample_node_parameters_zip_bart(tree, N, v, delta, phi,
                                  mu_residual, lambda_residual,
                                  alpha1, beta1, alpha2, beta2)
  
  return(list(tree = param$tree, mu = param$mu, lambda = param$lambda))
}


# 전체 트리에 대한 잠재 변수 갱신
update_latent_variables_from_all_trees <- function(tree_list, N, v,
                                                   alpha1, beta1, alpha2, beta2) {
  n <- length(N)
  delta_all <- numeric(n)
  phi_all <- numeric(n)
  
  for (tree in tree_list) {
    latent <- update_latent_variables(tree, N, v, alpha1, beta1, alpha2, beta2)
    delta_all <- delta_all + latent$delta
    phi_all <- phi_all + latent$phi
  }
  
  return(list(delta = delta_all / length(tree_list),
              phi = phi_all / length(tree_list)))
}

#mu, lambda sampling
sample_node_parameters_zip_bart <- function(tree, N, v, delta, phi,
                                   mu_residual, lambda_residual,
                                   alpha1, beta1, alpha2, beta2) {
  mu <- rep(NA, length(N))
  lambda <- rep(NA, length(N))
  node_id <- rep(NA, length(N))
  
  for (i in seq_along(tree$nodes)) {
    node <- tree$nodes[[i]]
    
    if (node$is_terminal) {
      indices <- node$data_indices
      
      # 🔧 잔차 기반 sufficient statistics
      sum_mu_resid <- sum(mu_residual[indices])
      sum_lambda_resid <- sum(lambda_residual[indices])
      sum_v <- sum(v[indices])  # 여전히 노출은 그대로 사용
      
      # mu 샘플링 (zero-inflation part)
      mu_sample <- rgamma(1,
                          shape = alpha1 + sum_mu_resid,
                          rate  = beta1 + length(indices))  # or another rate choice
      
      # lambda 샘플링 (Poisson part)
      lambda_sample <- rgamma(1,
                              shape = alpha2 + sum_lambda_resid,
                              rate  = beta2 + sum_v)
      
      # 노드에 저장
      node$mu <- mu_sample
      node$lambda <- lambda_sample
      
      # 트리에도 반영
      tree$nodes[[i]] <- node
      
      # 관측치별 벡터로 복사
      mu[indices] <- mu_sample
      lambda[indices] <- lambda_sample
      node_id[indices] <- node$node_id
    }
  }
  
  return(list(mu = mu, lambda = lambda, node_id = node_id, tree = tree))
}

# ZIP-BART 학습 함수 (m개 트리)
draw_zip_bart_model <- function(X, N, v, m = 50, n_iter = 100,
                                alpha1 = 1, beta1 = 5,
                                alpha2 = 1, beta2 = 2,
                                alpha_tree = 0.95, beta_tree = 1.0,
                                verbose = TRUE) {
  n <- length(N)
  tree_list <- vector("list", m)
  mu_matrix <- matrix(0, nrow = n, ncol = m)
  lambda_matrix <- matrix(0, nrow = n, ncol = m)
  
  delta <- rep(1, n)
  phi <- rexp(n, rate = 1 + 1)
  
  for (j in 1:m) {
    tree_list[[j]] <- create_initial_tree(X, N, v)
  }
  
  for (iter in 1:n_iter) {
    if (verbose) cat("[ZIP-BART Iteration", iter, "]\n")
    
    # Step 1: update latent variables
    latent <- update_latent_variables_from_all_trees(tree_list, N, v,
                                                     alpha1, beta1, alpha2, beta2)
    delta <- latent$delta
    phi <- latent$phi
    
    mu_target <- delta / phi
    lambda_target <- N / pmax(v, 1e-8)
    
    for (j in 1:m) {
      mu_minus_j <- rowSums(mu_matrix[, -j, drop = FALSE])
      lambda_minus_j <- rowSums(lambda_matrix[, -j, drop = FALSE])
      
      residual_mu <- mu_target - mu_minus_j
      residual_mu <- pmax(residual_mu, 1e-8)
      residual_mu[is.na(residual_mu) | is.infinite(residual_mu)] <- 1e-8
      residual_lambda <- lambda_target - lambda_minus_j
      residual_lambda <- pmax(residual_lambda, 1e-8)  # 음수/0 방지
      residual_lambda[is.na(residual_lambda) | is.infinite(residual_lambda)] <- 1e-8
      
      
      res <- update_zip_bart_tree(
        tree = tree_list[[j]], X = X, N = N, v = v,
        delta = delta, phi = phi,
        mu_residual = residual_mu, lambda_residual = residual_lambda,
        alpha1 = alpha1, beta1 = beta1,
        alpha2 = alpha2, beta2 = beta2,
        alpha_tree = alpha_tree, beta_tree = beta_tree
      )
      
      tree_list[[j]] <- res$tree
      mu_matrix[, j] <- res$mu
      lambda_matrix[, j] <- res$lambda
    }
  }
  
  return(list(tree_list = tree_list,
              mu_matrix = mu_matrix,
              lambda_matrix = lambda_matrix))
}

# ZIP-BART 예측 함수
predict_zip_bart <- function(tree_list, X_test, v_test) {
  m <- length(tree_list)
  n_test <- nrow(X_test)
  
  mu_mat <- matrix(0, nrow = n_test, ncol = m)
  lambda_mat <- matrix(0, nrow = n_test, ncol = m)
  
  for (j in 1:m) {
    tree <- tree_list[[j]]
    
    for (i in 1:n_test) {
      node <- find_terminal_node(tree, X_test[i, , drop = FALSE])
      
      # 안전하게 처리
      mu_ij <- if (!is.null(node$mu) && is.finite(node$mu)) node$mu else 1e-8
      lambda_ij <- if (!is.null(node$lambda) && is.finite(node$lambda)) node$lambda else 1e-8
      
      mu_mat[i, j] <- mu_ij
      lambda_mat[i, j] <- lambda_ij
    }
  }
  
  # ZIP 기대값 계산: 트리 합산
  zip_mean <- rowSums((mu_mat / (1 + mu_mat)) * lambda_mat) * v_test
  return(zip_mean)
}
