############################################ Tree Proposal Functions #######################################



tree_grow <- function(tree){
  
  # grow the tree
  # args: tree: current tree
  
  tree_index <- tree$index # index of tree terminal nodes for visits
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1] # index of tree termainal nodes for individuals
  tree_copy <- Clone(tree) # make a copy of the tree
  tree_grow_success <- TRUE # index of successfully grow the tree
  
  # randomly pick a terminal node for splitting
  terminal_node_name_all <- as.vector(tree$Get("name", filterFun = isLeaf))
  terminal_node_number <- length(terminal_node_name_all)
  sample_index <- sample(1:terminal_node_number, terminal_node_number, replace = FALSE)
  for (i in sample_index){
    terminal_node_name <- terminal_node_name_all[i]
    terminal_node <- FindNode(tree, paste(`terminal_node_name`))
    # if only the terminal node is available for splitting 
    # i.e., if the terminal node contains at least minimum number of individuals
    if (sum(tree_index_ind==terminal_node_name) > min_num) { 
      break 
    }
  }
  
  # randomly pick a split position 
  if (length(terminal_node$set) > 1){
    split_position <- sample(terminal_node$set, 1, replace = FALSE)
  }else { split_position <- terminal_node$set }
  terminal_node$split <- split_position
  
  # randomly assign the splitting rule 
  if (data_type[split_position] == "binary"){
    terminal_node$rule <- 0
    # binary covariate is unavailable for splitting later
    terminal_node$set <- terminal_node$set[terminal_node$set!=split_position]
  }else if (data_type[split_position] == "continuous"){
    # use uniform distribution to generate continuous splitting value
    epsilon <- 1e-10 # avoid numerical issue and empty terminal node
    available_value <- Xtree[tree_index==terminal_node_name,split_position]
    lower <- min(available_value) + epsilon
    upper <- max(available_value) - epsilon
    terminal_node$rule <- runif(1, lower, upper)
  }else if (data_type[split_position] == "categorical"){
    # randomly generate a non-empty subset of all available categories as splitting rule
    available_subset <- c()
    available_value <- sort(unique(Xtree[tree_index==terminal_node_name,split_position]))
    if (length(available_value) <= 1){
      tree_grow_success <- FALSE
    }else{
      for (c in 1:(length(available_value)-1)){
        # all subsets of available_value with c categories 
        subset <- combn(available_value, c) 
        available_subset <- c(available_subset, lapply(seq_len(ncol(subset)), function(i) subset[,i]))
      }
      index_rule <- sample(1:(2^(length(available_value))-2), 1, replace=F)
      terminal_node$rule <- available_subset[[index_rule]]
    }
  }else if (data_type[split_position] == "ordinal"){
    # randomly select an available level as tree splitting rule
    available_value <- sort(unique(Xtree[tree_index==terminal_node_name,split_position]))
    available_value <- available_value[1:(length(available_value)-1)] # avoid empty tree terminal node
    if (length(available_value) < 1){
      tree_grow_success <- FALSE
    }else if (length(available_value) == 1){
      terminal_node$rule <- available_value
    }else{ terminal_node$rule <- sample(available_value, 1, replace = FALSE) }
  }
  
  # add children to the split node
  terminal_node_info <- as.numeric(strsplit(as.character(terminal_node_name), "-")[[1]])
  terminal_node_depth <- terminal_node_info[1]; terminal_node_index <- terminal_node_info[2]
  children_node_depth <- terminal_node_depth + 1
  lchild_node_index <- 2*(terminal_node_index-1) + 1
  rchild_node_index <- 2*(terminal_node_index)
  
  # add left child 
  lchild_node_name <- paste(`children_node_depth`, "-", `lchild_node_index`, sep="")
  terminal_node$AddChild(paste(`lchild_node_name`))
  # add right child
  rchild_node_name <- paste(`children_node_depth`, "-", `rchild_node_index`, sep="")
  terminal_node$AddChild(paste(`rchild_node_name`))
  
  # update index of tree terminal nodes 
  data <- Xtree[tree_index==terminal_node_name,split_position]
  tree_index_subset <- tree_index[tree_index==terminal_node_name]
  if (data_type[split_position] == "binary"){
    tree_index_subset[which(data==terminal_node$rule)] <- lchild_node_name
    tree_index_subset[which(data!=terminal_node$rule)] <- rchild_node_name
  }else if (data_type[split_position] == "continuous"){
    tree_index_subset[which(data<=terminal_node$rule)] <- lchild_node_name
    tree_index_subset[which(data>terminal_node$rule)] <- rchild_node_name
  }else if (data_type[split_position] == "categorical"){
    tree_index_subset[which(data%in%terminal_node$rule)] <- lchild_node_name
    tree_index_subset[which(!(data%in%terminal_node$rule))] <- rchild_node_name
  }else if (data_type[split_position] == "ordinal"){
    tree_index_subset[which(data<=terminal_node$rule)] <- lchild_node_name
    tree_index_subset[which(data>terminal_node$rule)] <- rchild_node_name
  }
  tree_index[tree_index==terminal_node_name] <- tree_index_subset
  tree$index <- tree_index
  
  # tree will not sucessfully grow if exists small terminal node
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1]  
  if (sum(tree_index_ind==lchild_node_name) < min_num | sum(tree_index_ind==rchild_node_name) < min_num){
    tree_grow_success <- FALSE
  }
  
  # update the tree prior 
  terminal_node$node_prior <- 1-terminal_node$node_prior
  lchild_node <- FindNode(terminal_node, paste(`lchild_node_name`))
  rchild_node <- FindNode(terminal_node, paste(`rchild_node_name`))
  lchild_node$node_prior <- 1-a*(1+children_node_depth)^(-b)
  rchild_node$node_prior <- 1-a*(1+children_node_depth)^(-b)
  tree$prior <- prod(as.numeric(tree$Get("node_prior")))
  
  # update the available covariate split position at children nodes
  # note that binary covariate is unavailable for splitting later already
  lchild_node$set <- terminal_node$set
  rchild_node$set <- terminal_node$set
  if (data_type[split_position] %in% c("categorical","ordinal")){  
    # if the children node only contains one category or level then it is unavailable for splitting later 
    if (length(sort(unique(Xtree[tree_index==lchild_node_name,split_position]))) <=1 ){
      lchild_node$set <- lchild_node$set[lchild_node$set!=split_position] 
    }
    if (length(sort(unique(Xtree[tree_index==rchild_node_name,split_position]))) <=1 ){
      rchild_node$set <- rchild_node$set[rchild_node$set!=split_position] 
    }
  }
  
  # update the splitting node name 
  tree$old_node <- c(terminal_node_name)
  tree$new_node <- c(lchild_node_name, rchild_node_name)
  tree_copy$old_node <- c(terminal_node_name)
  tree_copy$new_node <- c(terminal_node_name)
  
  # update the tree if it was successfully grown
  if (tree_grow_success){
    tree$update_success <- tree_grow_success
    return(tree) 
  }else { 
    tree_copy$update_success <- tree_grow_success
    return(tree_copy) 
  }
}



tree_prune <- function(tree){
  
  # prune the tree
  # args: tree: current tree
  
  tree_index <- tree$index # index of tree terminal nodes 
  
  # randomly pick a parent node with two terminal children 
  terminal_node_name_all <- as.vector(tree$Get("name", filterFun = isLeaf))
  terminal_node_number <- length(terminal_node_name_all)
  sample_index <- sample(1:terminal_node_number, terminal_node_number, replace = FALSE)
  for (i in sample_index){
    terminal_node_name <- terminal_node_name_all[i]
    terminal_node <- FindNode(tree, paste(`terminal_node_name`))
    # if only the terminal node's sibling is also terminal node
    # i.e., its parent's children are both terminal nodes 
    if (names(terminal_node$siblings) %in% terminal_node_name_all){
      parent_node <- terminal_node$parent
      parent_node_name <- parent_node$name
      break
    }
  }
  
  # update index of tree terminal nodes 
  children_node_name <- names(parent_node$children)
  tree_index[tree_index %in% children_node_name] <- parent_node_name
  tree$index <- tree_index
  
  # collaspse the terminal nodes 
  lchild_node_name <- children_node_name[1]
  rchild_node_name <- children_node_name[2]
  parent_node$RemoveChild(lchild_node_name) 
  parent_node$RemoveChild(rchild_node_name) 
  # if the split position is binary covariate, then add it back 
  if (data_type[parent_node$split] == "binary"){
    parent_node$set <- sort(c(parent_node$split, parent_node$set))
  }
  parent_node$split <- c()
  parent_node$rule <- c()
  
  # update the tree prior
  parent_node$node_prior <- 1-parent_node$node_prior 
  tree$prior <- prod(as.numeric(tree$Get("node_prior")))
  
  # update the collapsing node name 
  tree$old_node <- c(lchild_node_name, rchild_node_name)
  tree$new_node <- c(parent_node_name)
  
  tree$update_success <- TRUE # tree prune always success
  return(tree)
}



tree_change_split <- function(tree){
  
  # change the split position of the tree
  # args: tree: current tree
  
  tree_index <- tree$index # index of tree terminal nodes for visits
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1] # index of tree termainal nodes for individuals
  tree_copy <- Clone(tree) # make a copy of the tree
  tree_change_success <- TRUE # index of successfully change the tree
  
  # randomly pick an internal node for changing the split
  internal_node_name <- as.vector(sample(tree$Get("name", filterFun = isNotLeaf), 1, replace = FALSE))
  internal_node <- FindNode(tree, paste(`internal_node_name`))
  terminal_node_name <- internal_node$Get("name", filterFun = isLeaf)
  
  # previous split position and available set 
  previous_split <- internal_node$split
  previous_set <- internal_node$set
  # all the splitting position for descendant internal nodes 
  descendant_split <- internal_node$Get("split", filterFun = isNotLeaf)
  descendant_split <- descendant_split[names(descendant_split)!=internal_node_name]
  
  # if descendant has binary split position, then remove it from available set
  # otherwise, the available set is the same as previous available set 
  if ("binary" %in% data_type[descendant_split]){
    binary_index <- which(data_type == "binary")
    available_set <- previous_set[!(previous_set%in%binary_index)]
  }else { available_set <- previous_set }
  # previous split position is not available for splitting 
  available_set <- available_set[available_set!=previous_split]
  
  # start to change the split
  if (length(available_set) == 0){
    tree_change_success <- FALSE
  }else {
    
    # randomly pick a split position 
    if (length(available_set) > 1){
      split_position <- sample(available_set, 1, replace = FALSE)
    }else { split_position <- available_set }
    # update the split position and available set 
    internal_node$split <- split_position 
    internal_node$set <- sort(unique(c(previous_split, previous_set)))
    # determine the lower and upper bound of the available set 
    available_set_list <- tree_available_set_top_down(internal_node, tree_index, split_position)
    lower <- available_set_list$lower
    upper <- available_set_list$upper 
    
    # randomly assign the splitting rule 
    if (data_type[split_position] == "binary"){
      # update the binary splitting rule 
      internal_node$rule <- 0
      # binary covariate is unavailable for splitting later
      internal_node$set <- internal_node$set[internal_node$set!=split_position]
    }else if (data_type[split_position] == "continuous"){
      # update the continuous splitting rule 
      if (lower > upper){
        tree_change_success <- FALSE
      }else { internal_node$rule <- runif(1, lower, upper) }
    }else if (data_type[split_position] == "categorical"){
      # update the categorical splitting rule 
      if (all(upper %in% lower) & !(identical(lower, upper))){
        tree_change_success <- FALSE
      }else{
        available_subset <- c()
        available_value <- sort(unique(Xtree[tree_index%in%terminal_node_name,split_position]))
        if (length(available_value) <= 1){
          tree_change_success <- FALSE
        }else {
          # all the available subsets without considering lower and upper bound
          for (c in 1:(length(available_value)-1)){
            subset <- combn(available_value, c)
            available_subset <- c(available_subset, lapply(seq_len(ncol(subset)), function(i) subset[,i]))
          }
          # delete the subsets outside the lower and upper bound
          available_index <- c()
          for (index in 1:(2^(length(available_value))-2)){
            if (all(lower%in%available_subset[[index]]) & all(available_subset[[index]]%in%upper)){
              available_index <- c(available_index, index)
            } 
          }
          # randomly select an available subset 
          if (length(available_index) < 1){
            tree_change_success <- FALSE
          }else if (length(available_index) == 1){
            index_rule <- available_index
            internal_node$rule <- available_subset[[index_rule]]
          }else if (length(available_index) > 1){
            index_rule <- sample(available_index, 1, replace = FALSE)
            internal_node$rule <- available_subset[[index_rule]]
          }
        }
      }
    }else if (data_type[split_position] == "ordinal"){
      # update the ordinal splitting rule 
      if (lower > upper){
        tree_change_success <- FALSE
      }else { 
        if (lower == upper){
          internal_node$rule <- lower
        }else { 
          internal_node$rule <- sample(lower:upper, 1, replace = FALSE)
        }
      }
    }
    
    # update index of tree terminal nodes 
    tree_index_update_list <- tree_index_update(internal_node, tree_index, tree_change_success)
    tree_index <- tree_index_update_list$tree_index
    tree_change_success <- tree_index_update_list$tree_change_success
    tree$index <- tree_index
    
    # update the available set of descendant nodes 
    tree_descendant_set_update(internal_node, tree_index, split_position, previous_split)
  }
  
  # tree will not sucessfully change if exists small terminal node
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1]
  terminal_node_name <- tree$Get("name", filterFun = isLeaf)
  for (node_name in terminal_node_name){
    if (sum(tree_index_ind == node_name) < min_num){
      tree_change_success <- FALSE
    }
  }
  
  # update the tree if it was successfully changed
  if (tree_change_success){
    tree$update_success <- tree_change_success
    return(tree) 
  }else { 
    tree_copy$update_success <- tree_change_success
    return(tree_copy) 
  }
}  



tree_change_rule <- function(tree){
  
  # change the splitting rule of the tree
  # args: tree: current tree
  
  tree_index <- tree$index # index of tree terminal nodes for visits
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1] # index of tree termainal nodes for individuals
  tree_copy <- Clone(tree) # make a copy of the tree
  tree_change_success <- TRUE # index of successfully change the tree
  
  # randomly pick a split position for changing the rule
  split_position_all <- tree$Get("split", filterFun = isNotLeaf)
  split_position_all <- split_position_all[which(data_type[split_position_all]!="binary")]
  if (length(split_position_all) < 1){
    tree_change_success <- FALSE
    split_position <- tree$root$split
  }else if (length(split_position_all) == 1){
    split_position <- split_position_all
  }else { split_position <- sample(split_position_all, 1, replace = FALSE) }
  internal_node_name_all <- sort(unique(names(split_position_all[split_position_all==split_position])))
  
  # update the splitting rule by the topological order of the tree internal nodes
  for (internal_node_name in internal_node_name_all){
    
    internal_node <- FindNode(tree, paste(`internal_node_name`))
    # path from root node to internal node
    node_along_path <- internal_node$path 
    node_along_path <- node_along_path[node_along_path!=internal_node_name]
    node_along_path <- intersect(node_along_path, internal_node_name_all)
    # determine the lower and upper bound of the available set 
    available_set_list <- tree_available_set_bottom_up(tree, internal_node, node_along_path, tree_index, split_position)
    lower <- available_set_list$lower
    upper <- available_set_list$upper 
    
    # update the splitting rule for each internal nodes 
    if (data_type[split_position] == "continuous"){
      # update the continuous splitting rule 
      if (lower > upper){
        tree_change_success <- FALSE
      }else { 
        internal_node$rule <- runif(1, lower, upper) 
      }
    }else if (data_type[split_position] == "categorical"){
      # update the categorical splitting rule 
      if (all(upper %in% lower) & !(identical(lower, upper))){
        tree_change_success <- FALSE
      }else{
        available_subset <- c()
        available_value <- sort(unique(Xtree[data_index,split_position]))
        # all the available subsets without considering lower and upper bound
        for (c in 1:(length(available_value)-1)){
          subset <- combn(available_value, c)
          available_subset <- c(available_subset, lapply(seq_len(ncol(subset)), function(i) subset[,i]))
        }
        # delete the subsets outside the lower and upper bound
        available_index <- c()
        for (index in 1:(2^(length(available_value))-2)){
          if (all(lower%in%available_subset[[index]]) & all(available_subset[[index]]%in%upper)){
            available_index <- c(available_index, index)
          } 
        }
        # randomly select an available subset 
        if (length(available_index) > 1){
          index_rule <- sample(available_index, 1, replace = FALSE)
        }else { index_rule <- available_index }
        internal_node$rule <- available_subset[[index_rule]]
      }
    }else if (data_type[split_position] == "ordinal"){
      # update the ordinal splitting rule 
      if (lower > upper){
        tree_change_success <- FALSE
      }else { 
        if (lower == upper){
          internal_node$rule <- lower
        }else { 
          internal_node$rule <- sample(lower:upper, 1, replace = FALSE)
        }
      }
    }
  }
  
  # update index of tree terminal nodes 
  tree_index_update_list <- tree_index_update(tree$root, tree_index, tree_change_success)
  tree_index <- tree_index_update_list$tree_index
  tree_change_success <- tree_index_update_list$tree_change_success
  tree$index <- tree_index
  
  # update the available set of descendant nodes 
  tree_descendant_set_update(tree$root, tree_index, split_position, split_position)
  
  # tree will not sucessfully change if exists small terminal node
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1]
  terminal_node_name <- tree$Get("name", filterFun = isLeaf)
  for (node_name in terminal_node_name){
    if (sum(tree_index_ind == node_name) < min_num){
      tree_change_success <- FALSE
    }
  }
  
  # update the tree if it was successfully changed
  if (tree_change_success){
    tree$update_success <- tree_change_success
    return(tree) 
  }else { 
    tree_copy$update_success <- tree_change_success
    return(tree_copy) 
  }
}



tree_swap <- function(tree){
  
  # swap the tree
  # args: tree: current tree
  
  tree_index <- tree$index # index of tree terminal nodes for visits
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1] # index of tree termainal nodes for individuals
  tree_copy <- Clone(tree) # make a copy of the tree
  tree_change_success <- TRUE # index of successfully change the tree
  
  # randomly pick a pair of parent-child internal nodes
  internal_node_name_all <- as.vector(tree$Get("name", filterFun = isNotLeaf))
  internal_node_name_all <- internal_node_name_all[internal_node_name_all!=tree$root$name]
  internal_node_number <- length(internal_node_name_all)
  sample_index <- sample(1:internal_node_number, internal_node_number, replace = FALSE)
  for (i in sample_index){
    child_node_name <- internal_node_name_all[i]
    child_node <- FindNode(tree, paste(`child_node_name`))
    parent_node <- child_node$parent
    parent_node_name <- parent_node$name
    # if only their split positions are different 
    if (parent_node$split != child_node$split){
      break
    }
  }
  
  # previous split positions and rules
  parent_split <- parent_node$split
  parent_rule <- parent_node$rule
  child_split <- child_node$split
  child_rule <- child_node$rule
  sibling_node <- child_node$siblings
  sibling_split <- sibling_node$split
  sibling_rule <- sibling_node$rule
  # to simplify the following if condition,
  # if the sibling node is a terminal node,
  # then assign it zero split position and rule 
  if (is.null(sibling_split)){
    sibling_split <- 0; sibling_rule <- 0
  }
  
  # swap the split position and rule 
  if (child_split == sibling_split & identical(child_rule, sibling_rule)){
    # if the child node and its sibling node share the same split position and rule
    parent_node$split <- child_split
    parent_node$rule <- child_rule
    child_node$split <- parent_split
    child_node$rule <- parent_rule
    sibling_node$split <- parent_split
    sibling_node$rule <- parent_rule
  }else {
    parent_node$split <- child_split
    parent_node$rule <- child_rule
    child_node$split <- parent_split
    child_node$rule <- parent_rule
  }
  
  # update index of tree terminal nodes 
  tree_index_update_list <- tree_index_update(parent_node, tree_index, tree_change_success)
  tree_index <- tree_index_update_list$tree_index
  tree_change_success <- tree_index_update_list$tree_change_success
  tree$index <- tree_index
  
  # check whether the tree can be successfully changed by 
  # compare the split rule with the lower and upper bound of the available value
  if (child_split == sibling_split & identical(child_rule, sibling_rule)){
    
    # update the available sets 
    parent_node$set <- sort(unique(c(parent_split, parent_node$set)))
    child_node$set <- sort(unique(c(child_split, child_node$set)))
    sibling_node$set <- sort(unique(c(sibling_split, sibling_node$set)))
    
    # determine the lower and upper bound of the available set for parent node
    parent_available_set_list <- tree_available_set_top_down(parent_node, tree_index, child_split)
    parent_lower <- parent_available_set_list$lower
    parent_upper <- parent_available_set_list$upper 
    if (data_type[child_split] == "binary"){
      parent_node$set <- parent_node$set[parent_node$set!=child_split]
    }else if (data_type[child_split] == "continuous"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (child_rule < parent_lower | child_rule > parent_upper){
        tree_change_success <- FALSE
      }
    }else if (data_type[child_split] == "categorical"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (all(child_rule %in% parent_lower) | all(parent_upper %in% child_rule)){
        tree_change_success <- FALSE
      }
    }else if (data_type[child_split] == "ordinal"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (child_rule < parent_lower | child_rule > parent_upper){
        tree_change_success <- FALSE
      }
    }
    
    # determine the lower and upper bound of the available set for child node
    child_available_set_list <- tree_available_set_top_down(child_node, tree_index, parent_split)
    child_lower <- child_available_set_list$lower
    child_upper <- child_available_set_list$upper 
    if (data_type[parent_split] == "binary"){
      child_node$set <- child_node$set[child_node$set!=parent_split]
    }else if (data_type[parent_split] == "continuous"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < child_lower | parent_rule > child_upper){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "categorical"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (all(parent_rule %in% child_lower) | all(child_upper %in% parent_rule)){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "ordinal"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < child_lower | parent_rule > child_upper){
        tree_change_success <- FALSE
      }
    }
    
    # determine the lower and upper bound of the available set for sibling node
    sibling_available_set_list <- tree_available_set_top_down(sibling_node, tree_index, parent_split)
    sibling_lower <- sibling_available_set_list$lower
    sibling_upper <- sibling_available_set_list$upper 
    if (data_type[parent_split] == "binary"){
      sibling_node$set <- sibling_node$set[sibling_node$set!=parent_split]
    }else if (data_type[parent_split] == "continuous"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < sibling_lower | parent_rule > sibling_upper){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "categorical"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (all(parent_rule %in% sibling_lower) | all(child_upper %in% sibling_rule)){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "ordinal"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < sibling_lower | parent_rule > sibling_upper){
        tree_change_success <- FALSE
      }
    }
  }else {
    
    # update the available sets 
    parent_node$set <- sort(unique(c(parent_split, parent_node$set)))
    child_node$set <- sort(unique(c(child_split, child_node$set)))
    
    # determine the lower and upper bound of the available set for parent node
    parent_available_set_list <- tree_available_set_top_down(parent_node, tree_index, child_split)
    parent_lower <- parent_available_set_list$lower
    parent_upper <- parent_available_set_list$upper 
    if (data_type[child_split] == "binary"){
      # if the child split at binary variable, 
      # and its sibling's descendant splits at the same position, 
      # then the tree is not able to be swapped 
      if (isLeaf(sibling_node)){
        sibling_descendant_split <- c()
      }else {
        sibling_descendant_split <- sibling_node$Get("split", filterFun = isNotLeaf)
      }
      if (child_split %in% sibling_descendant_split){
        tree_change_success <- FALSE
      }else{
        parent_node$set <- parent_node$set[parent_node$set!=child_split]
      }
    }else if (data_type[child_split] == "continuous"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (child_rule < parent_lower | child_rule > parent_upper){
        tree_change_success <- FALSE
      }
    }else if (data_type[child_split] == "categorical"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (all(child_rule %in% parent_lower) | all(parent_upper %in% child_rule)){
        tree_change_success <- FALSE
      }
    }else if (data_type[child_split] == "ordinal"){
      # if the child rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (child_rule < parent_lower | child_rule > parent_upper){
        tree_change_success <- FALSE
      }
    }
    
    # determine the lower and upper bound of the available set for child node
    child_available_set_list <- tree_available_set_top_down(child_node, tree_index, parent_split)
    child_lower <- child_available_set_list$lower
    child_upper <- child_available_set_list$upper 
    if (data_type[parent_split] == "binary"){
      child_node$set <- child_node$set[child_node$set!=parent_split]
    }else if (data_type[parent_split] == "continuous"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < child_lower | parent_rule > child_upper){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "categorical"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (all(parent_rule %in% child_lower) | all(child_upper %in% parent_rule)){
        tree_change_success <- FALSE
      }
    }else if (data_type[parent_split] == "ordinal"){
      # if the parent rule is outside the lower and upper bound 
      # then the tree is not able to be swapped 
      if (parent_rule < child_lower | parent_rule > child_upper){
        tree_change_success <- FALSE
      }
    }
  }
  
  # update the available set of descendant nodes 
  if (child_split == sibling_split & identical(child_rule, sibling_rule)){
    tree_descendant_set_update(parent_node, tree_index, child_split, parent_split)
    tree_descendant_set_update(child_node, tree_index, parent_split, child_split)
    tree_descendant_set_update(sibling_node, tree_index, parent_split, sibling_split)
  }else{
    tree_descendant_set_update(parent_node, tree_index, child_split, parent_split)
    tree_descendant_set_update(child_node, tree_index, parent_split, child_split)
  }
  
  # tree will not sucessfully change if exists small terminal node
  tree_index_ind <- matrix(tree_index, nrow=n, ncol=max(J))[,1]
  terminal_node_name <- tree$Get("name", filterFun = isLeaf)
  for (node_name in terminal_node_name){
    if (sum(tree_index_ind == node_name) < min_num){
      tree_change_success <- FALSE
    }
  }
  
  # update the tree if it was successfully changed
  if (tree_change_success){
    tree$update_success <- tree_change_success
    return(tree) 
  }else { 
    tree_copy$update_success <- tree_change_success
    return(tree_copy) 
  }
}



tree_available_set_top_down <- function(internal_node, tree_index, split_position){
  
  # determine the available set of the tree top down (for tree_change_split and tree_swap)
  # args: internal_node: the internal node
  #       tree_index: index of current tree terminal node 
  #       split_position: the split position of the internal node 
  
  internal_node_name <- internal_node$name
  terminal_node_name <- internal_node$Get("name", filterFun = isLeaf)
  lower <- NA; upper <- NA # initialize lower and upper bound with NA values
  
  # the splitting position and rule at children node branch 
  children_node_name <- names(internal_node$children) 
  lchild_node_name <- children_node_name[1]
  lchild_node <- FindNode(internal_node, paste(`lchild_node_name`))
  lchild_node_descendant_split <- lchild_node$Get("split", filterFun = isNotLeaf)
  lchild_node_descendant_rule <- lchild_node$Get("rule", filterFun = isNotLeaf)
  rchild_node_name <- children_node_name[2]
  rchild_node <- FindNode(internal_node, paste(`rchild_node_name`))
  rchild_node_descendant_split <- rchild_node$Get("split", filterFun = isNotLeaf)
  rchild_node_descendant_rule <- rchild_node$Get("rule", filterFun = isNotLeaf)
  
  if (data_type[split_position] == "continuous"){
    epsilon <- 1e-10 # avoid numerical issue and empty terminal node
    available_value <- Xtree[tree_index%in%terminal_node_name,split_position]
    lower <- min(available_value) + epsilon
    upper <- max(available_value) - epsilon
    # update the lower bound by the left descendant node splitting rule
    if (split_position %in% lchild_node_descendant_split){
      split_position_index <- which(lchild_node_descendant_split == split_position)
      for (i in split_position_index){ 
        lower <- max(lower, lchild_node_descendant_rule[[i]]) 
      }
    }
    # update the upper bound by the right descendant node splitting rule
    if (split_position %in% rchild_node_descendant_split){
      split_position_index <- which(rchild_node_descendant_split == split_position)
      for (i in split_position_index){ 
        upper <- min(upper, rchild_node_descendant_rule[[i]]) 
      }
    }
  }else if (data_type[split_position] == "categorical"){
    # randomly generate a non-empty subset of all available categories as splitting rule
    available_value <- sort(unique(Xtree[tree_index%in%terminal_node_name,split_position]))
    lower <- c()
    upper <- available_value
    # update the lower bound by the left descendant node splitting rule
    if (split_position %in% lchild_node_descendant_split){
      lchild_node_descendant_name <- lchild_node$Get("name", filterFun = isLeaf)
      lchild_node_descendant_value <- sort(unique(Xtree[tree_index%in%lchild_node_descendant_name,split_position]))
      lower <- union(lower, lchild_node_descendant_value)
    }
    # update the upper bound by the right descendant node splitting rule
    if (split_position %in% rchild_node_descendant_split){
      rchild_node_descendant_name <- rchild_node$Get("name", filterFun = isLeaf)
      rchild_node_descendant_value <- sort(unique(Xtree[tree_index%in%rchild_node_descendant_name,split_position]))
      upper <- setdiff(upper, rchild_node_descendant_value)
    }
  }else if (data_type[split_position] == "ordinal"){
    available_value <- Xtree[tree_index%in%terminal_node_name,split_position]
    lower <- min(available_value)
    upper <- max(available_value)
    # update the lower bound by the left descendant node splitting rule
    if (split_position %in% lchild_node_descendant_split){
      split_position_index <- which(lchild_node_descendant_split == split_position)
      for (i in split_position_index){ 
        lower <- max(lower, lchild_node_descendant_rule[[i]]) 
      }
    }
    # update the upper bound by the right descendant node splitting rule
    if (split_position %in% rchild_node_descendant_split){
      split_position_index <- which(rchild_node_descendant_split == split_position)
      for (i in split_position_index){ 
        upper <- min(upper, rchild_node_descendant_rule[[i]]) 
      }
    }
  }
  
  # return the lower and upper bound of the available set
  available_set_list <- list(lower=lower, upper=upper)
  return(available_set_list)
}



tree_available_set_bottom_up <- function(tree, internal_node, node_along_path, tree_index, split_position){
  
  # determine the available set of the tree bottom up (for tree_change_rule)
  # args: tree: current tree
  #       internal_node: the internal node
  #       node_along_path: the related node along the path from internal node to root node 
  #       tree_index: index of current tree terminal node 
  #       split_position: the split position of the internal node 
  
  internal_node_name = internal_node$name
  path_length <- length(node_along_path) # number of internal nodes along the path 
  
  if (data_type[split_position] == "continuous"){
    epsilon <- 1e-10 # avoid numerical issue and empty terminal node
    available_value <- Xtree[data_index,split_position]
    lower <- min(available_value) + epsilon
    upper <- max(available_value) - epsilon
    if (path_length >= 1){
      # if there exists ancestor node share the same split position with the internal node 
      for (node_name in node_along_path){
        node <- FindNode(tree, paste(`node_name`))
        children_node_name <- names(node$children) 
        lchild_node_name <- children_node_name[1]
        # lchild_node <- FindNode(node, paste(`lchild_node_name`))
        lchild_node <- FindNode(tree, paste(`lchild_node_name`))
        if (internal_node_name %in% lchild_node$Get("name")){
          # if the internal node is in the left child branch of the node 
          upper <- min(upper, node$rule)
        }else {
          # if the internal node is in the right child branch of the node 
          lower <- max(lower, node$rule)
        }
      }
    }else {
      # if there does not exist ancestor node share the same split position with the internal node 
      step <- 0.05 # step for uniform random walk
      previous_rule <- internal_node$rule
      lower <- max(lower, previous_rule-step) 
      upper <- min(upper, previous_rule+step)
    }
  }else if (data_type[split_position] == "categorical"){
    available_value <- sort(unique(Xtree[data_index,split_position]))
    lower <- c()
    upper <- available_value
    if (path_length >= 1){
      # if there exists ancestor node share the same split position with the internal node 
      # otherwise, don't need to update the lower and upper bound of the categorical covariate
      for (node_name in node_along_path){
        node <- FindNode(tree, paste(`node_name`))
        children_node_name <- names(node$children) 
        lchild_node_name <- children_node_name[1]
        lchild_node <- FindNode(tree, paste(`lchild_node_name`))
        if (internal_node_name %in% lchild_node$Get("name")){
          # if the internal node is in the left child branch of the node 
          upper <- intersect(upper, node$rule)
        }else {
          # if the internal node is in the right child branch of the node 
          upper <- setdiff(upper, node$rule)
        }
      }
    }
  }else if (data_type[split_position] == "ordinal"){
    available_value <-  sort(unique(Xtree[data_index,split_position]))
    lower <- min(available_value)
    upper <- max(available_value)
    if (path_length >= 1){
      # if there exists ancestor node share the same split position with the internal node 
      # otherwise, don't need to update the lower and upper bound of the ordinal covariate
      for (node_name in node_along_path){
        node <- FindNode(tree, paste(`node_name`))
        children_node_name <- names(node$children) 
        lchild_node_name <- children_node_name[1]
        lchild_node <- FindNode(tree, paste(`lchild_node_name`))
        if (internal_node_name %in% lchild_node$Get("name")){
          # if the internal node is in the left child branch of the node 
          upper <- min(upper, node$rule)
        }else {
          # if the internal node is in the right child branch of the node 
          lower <- max(lower, node$rule)
        }
      }
    }
  }
  
  # return the lower and upper bound of the available set
  available_set_list <- list(lower=lower, upper=upper)
  return(available_set_list)
}



tree_index_update <- function(internal_node, tree_index, tree_change_success){
  
  # update the index of tree terminal nodes
  # args: internal_node: the internal node
  #       tree_index: index of current tree terminal node 
  #       tree_change_success: whether the tree can be successfully changed 
  
  subtree <- Clone(internal_node) 
  internal_node_name <- subtree$name
  terminal_node_name_all <- subtree$Get("name", filterFun = isLeaf)
  data <- Xtree[tree_index%in%terminal_node_name_all,]
  tree_index_subset <- tree_index[tree_index%in%terminal_node_name_all]
  
  # update the subtree index for each terminal nodes 
  for (terminal_node_name in terminal_node_name_all){
    
    # index of data in the terminal node (TRUE/FALSE)
    terminal_node <- FindNode(subtree, paste(`terminal_node_name`))
    terminal_node_data <- rep(TRUE, length(tree_index_subset)) 
    # path from internal node to terminal node 
    node_along_path <- terminal_node$path 
    node_along_path <- node_along_path[node_along_path!=internal_node_name]
    path_length <- length(node_along_path)
    # split position and rule along the path 
    split_along_path <- rev(terminal_node$Get("split", traversal = "ancestor"))
    rule_along_path <- rev(terminal_node$Get("rule", traversal = "ancestor"))
    split_along_path <- split_along_path[!is.na(split_along_path)]
    rule_along_path <- rule_along_path[!is.na(rule_along_path)]
    # split direction (left or right) along the path
    direction_along_path <- c()
    for (node_name in node_along_path){
      node_index <- as.numeric(strsplit(as.character(node_name), "-")[[1]])[2]
      if (node_index %% 2 != 0){
        # if it is the left child of the previous node along the path 
        direction_along_path <- c(direction_along_path, "left")
      }else { 
        # if it is the right child of the previous node along the path 
        direction_along_path <- c(direction_along_path, "right") 
      }
    }
    
    # update the index of data in the terminal node along the path 
    # i.e., a data point is in the terminal node if it satisfies all the rules along the path
    for (i in 1:path_length){
      if (direction_along_path[i] == "left"){
        if (data_type[split_along_path[i]] == "binary"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] == rule_along_path[[i]])
        }else if (data_type[split_along_path[i]] == "continuous"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] <= rule_along_path[[i]])
        }else if (data_type[split_along_path[i]] == "categorical"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] %in% rule_along_path[[i]])
        }else if (data_type[split_along_path[i]] == "ordinal"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] <= rule_along_path[[i]])
        }
      }else if (direction_along_path[i] == "right"){
        if (data_type[split_along_path[i]] == "binary"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] != rule_along_path[[i]])
        }else if (data_type[split_along_path[i]] == "continuous"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] > rule_along_path[[i]])
        }else if (data_type[split_along_path[i]] == "categorical"){
          terminal_node_data <- terminal_node_data & (!(data[,split_along_path[i]] %in% rule_along_path[[i]]))
        }else if (data_type[split_along_path[i]] == "ordinal"){
          terminal_node_data <- terminal_node_data & (data[,split_along_path[i]] > rule_along_path[[i]])
        }
      }
    }
    
    # if tree change will result in empty terminal node,
    # then the tree will not be successfully changed
    if (sum(terminal_node_data) < 1){
      tree_change_success <- FALSE
    }else{
      # update the index of the current tree terminal node
      tree_index_subset[terminal_node_data] <- terminal_node_name
    }
  }
  
  # update the index of all tree terminal nodes 
  tree_index[tree_index%in%terminal_node_name_all] <- tree_index_subset
  # return the updated tree index 
  tree_index_update_list <- list(tree_index=tree_index, tree_change_success=tree_change_success)
  return(tree_index_update_list)
}



tree_descendant_set_update <- function(internal_node, tree_index, split_position, previous_split){
  
  # update the available set for the descendant of the internal node 
  # args: internal_node: the internal node
  #       tree_index: index of current tree terminal node 
  #       split_position: the split position of the internal node 
  #       previous_split: the previous split position of the internal node 
  
  internal_node_name <- internal_node$name
  descendant_node_name <- internal_node$Get("name")
  
  if (data_type[split_position] %in% c("binary", "categorical", "ordinal")){
    # binary, categorical or ordinal covariate is unavailable for splitting later
    # if the descendant node contains only one available category
    for (node_name in descendant_node_name){
      descendant_node <- FindNode(internal_node, paste(`node_name`))
      descendant_terminal_node_name <- descendant_node$Get("name", filterFun=isLeaf)
      if (length(sort(unique(Xtree[tree_index%in%descendant_terminal_node_name,split_position]))) <=1){
        descendant_node$set <- descendant_node$set[descendant_node$set!=split_position]
      }
    }
  }
  
  if (data_type[previous_split] %in% c("binary", "categorical", "ordinal")){
    # if previous splitting position is binary, categorical or ordinal covariate
    # then add it back to the available set of descendant nodes 
    # which contains more than one available categories
    for (node_name in descendant_node_name){
      descendant_node <- FindNode(internal_node, paste(`node_name`))
      descendant_terminal_node_name <- descendant_node$Get("name", filterFun=isLeaf)
      if (length(sort(unique(Xtree[tree_index%in%descendant_terminal_node_name,previous_split]))) > 1){
        descendant_node$set <- sort(unique(c(previous_split, descendant_node$set)))
      }
    }
  }
}



############################################ Tree Update Functions #######################################



update_tree <- function(tree, sigma2_omega, rho, sigma2, etheta, Btheta){
  
  # update tree
  # args: tree: current tree
  #       sigma2_omega, rho, sigma2: parameters
  #       etheta, Btheta: hyper-parameters
  # returns: tree: updated tree
  
  ratio <- -Inf
  
  # propose the new candidate tree
  propose_tree_list <- propose_tree(Clone(tree))
  tree_new <- propose_tree_list$tree_new
  update_step <- propose_tree_list$update_step
  update_success <- propose_tree_list$update_success
  
  # propose the new tree-specific parameters based on training data
  if (update_success){
    propose_pars_list <- propose_pars(update_step, Clone(tree), Clone(tree_new), sigma2_omega, rho, sigma2, etheta, Btheta)
    sigma2_omega_new <- propose_pars_list$sigma2_omega
    rho_new <- propose_pars_list$rho
    
    # tree prior and proposal 
    tree_prior <- tree$prior
    tree_new_prior <- tree_new$prior
    tree_proposal <- propose_tree_list$tree_proposal
    tree_new_proposal <- propose_tree_list$tree_new_proposal
    
    # tree log-likelihood based on testing data 
    # tree_logll <- logll_tree(Clone(tree), sigma2_omega, rho, sigma2, etheta, Btheta, frac=1-xi)
    # tree_new_logll <- logll_tree(Clone(tree_new), sigma2_omega_new, rho_new, sigma2, etheta, Btheta, frac=1-xi)
    R_tilde <- tree$index; L <- length(unique(R_tilde))-1
    R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1]
    R_index <- sort(unique(R)) # unique index of tree terminal nodes
    R_num <- as.numeric(factor(R, levels=R_index)) # numerical index of the tree terminal node 
    tree_logll <- logll_tree_rcpp(R_num, sigma2_omega, rho, sigma2, etheta, Btheta, L, S_tilde, J, Kappa, t, bX, y, frac=1-xi)
    
    R_tilde_new <- tree_new$index; L_new  <- length(unique(R_tilde_new))-1
    R_new <- matrix(R_tilde_new, nrow=n, ncol=max(J))[,1]
    R_new_index <- sort(unique(R_new)) # unique index of tree terminal nodes
    R_new_num <- as.numeric(factor(R_new, levels=R_new_index)) # numerical index of the tree terminal node 
    tree_new_logll <- logll_tree_rcpp(R_new_num, sigma2_omega_new, rho_new, sigma2, etheta, Btheta, L_new, S_tilde, J, Kappa, t, bX, y, frac=1-xi)
    
    # acceptance ratio of the new candidate tree 
    ratio <- tree_new_logll - tree_logll + log(tree_new_prior) - log(tree_prior) + log(tree_proposal) - log(tree_new_proposal) 
  }
  
  # tree marginal log-posterior 
  if (update_success){
    # tree_logpost <- logpost_tree(Clone(tree), sigma2_omega, rho, sigma2, etheta, Btheta)
    # tree_new_logpost <- logpost_tree(Clone(tree_new), sigma2_omega_new, rho_new, sigma2, etheta, Btheta)
    tree_logpost <- logpost_tree_rcpp(Clone(tree), sigma2_omega, rho, sigma2, etheta, Btheta)
    tree_new_logpost <- logpost_tree_rcpp(Clone(tree_new), sigma2_omega_new, rho_new, sigma2, etheta, Btheta)
  }else{
    # tree_logpost <- logpost_tree(Clone(tree), sigma2_omega, rho, sigma2, etheta, Btheta)
    tree_logpost <- logpost_tree_rcpp(Clone(tree), sigma2_omega, rho, sigma2, etheta, Btheta)
  }
  
  # accept or reject
  if (log(runif(1)) < ratio & update_success){ 
    tree <- tree_new
    tree$accept <- TRUE; print("accept!")
    tree$update_step <- update_step 
    tree$logpost <- tree_new_logpost
    tree$sigma2_omega <- sigma2_omega_new
    tree$rho <- rho_new
  }else { 
    tree$accept <- FALSE; print("reject!")
    tree$update_step <- update_step
    tree$logpost <- tree_logpost
    tree$sigma2_omega <- sigma2_omega
    tree$rho <- rho
  } 
  
  return(tree)
}



propose_tree <- function(tree){
  
  # propose the new candidate tree 
  # args: tree: current tree
  # returns: propose_tree_list:
  #          tree_new: the proposed new tree
  #          update_step: the step of update
  #          update_success: indicator whether update successfully
  #          tree_proposal: the probability of propose the tree
  #          tree_new_proposa: the probability of propose the new tree
  
  steps <- c("grow", "prune", "change split", "change rule", "swap") # name of possible steps
  num_steps <- length(steps) # number of possible steps
  depth <- tree$height-1 # depth of the tree
  
  # propose the new candidate tree
  if (depth == 0){
    
    # case 1: depth == 0, tree can only grow
    update_step <- c("grow") 
    tree_new <- tree_grow(Clone(tree)); print("grow!")
    tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
    tree_new_proposal <- 1
    if (tree_new_split_num > 1){
      tree_proposal <- 1/num_steps
    }else { tree_proposal <- 1/(num_steps-1) }
    
  }else if (depth >= depth_max){
    
    # case 2: depth >= depth_max, tree can not grow      
    steps <- steps[steps != "grow"] 
    tree_split_num <- length(unique(tree$Get("split", filterFun = isNotLeaf)))
    if (tree_split_num > 1){
      # if there are at least two different split positions
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!")  
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1) * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!")
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "swap"){
        tree_new <- tree_swap(Clone(tree)); print("swap!")
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }
    }else {
      # if there is only one split position
      steps <- steps[steps != "swap"] 
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!")  
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-2) * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!")
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }
    }
    
  }else if (depth == (depth_max-1)){
    
    # case 3: depth = (depth_max-1), new tree can not grow
    tree_split_num <- length(unique(tree$Get("split", filterFun = isNotLeaf)))
    if (tree_split_num > 1){
      # if there are at least two different split positions
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "grow"){ 
        tree_new <- tree_grow(Clone(tree)); print("grow!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps * 1/length(tree$Get("name", filterFun = isLeaf))
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2)
        }else { tree_proposal <- 1/(num_steps-2) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2) }   
      }else if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-2) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-2) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-2) }
      }else if (update_step == "swap"){
        tree_new <- tree_swap(Clone(tree)); print("swap!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-2) }
      }
    }else {
      # if there is only one split position
      steps <- steps[steps != "swap"] 
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "grow"){ 
        tree_new <- tree_grow(Clone(tree)); print("grow!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1) * 1/length(tree$Get("name", filterFun = isLeaf))
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2)
        }else { tree_proposal <- 1/(num_steps-2) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2) }   
      }else if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1) * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-2) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-2) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-2) }
      }
    }
    
  }else {
    
    # case 4: otherwise, no constraint on the tree proposal 
    tree_split_num <- length(unique(tree$Get("split", filterFun = isNotLeaf)))
    if (tree_split_num > 1){
      # if there are at least two different split positions
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "grow"){ 
        tree_new <- tree_grow(Clone(tree)); print("grow!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps * 1/length(tree$Get("name", filterFun = isLeaf))
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2)
        }else { tree_proposal <- 1/(num_steps-1) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2) }   
      }else if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "swap"){
        tree_new <- tree_swap(Clone(tree)); print("swap!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/num_steps
        if (tree_new_split_num > 1){
          tree_proposal <- 1/num_steps
        }else { tree_proposal <- 1/(num_steps-1) }
      }
    }else {
      # if there is only one split position
      steps <- steps[steps != "swap"] 
      update_step <- sample(steps, 1, replace = FALSE) 
      if (update_step == "grow"){ 
        tree_new <- tree_grow(Clone(tree)); print("grow!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1) * 1/length(tree$Get("name", filterFun = isLeaf))
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2)
        }else { tree_proposal <- 1/(num_steps-1) * 1/floor(length(tree_new$Get("name", filterFun = isLeaf))/2) }   
      }else if (update_step == "prune"){ 
        tree_new <- tree_prune(Clone(tree)); print("prune!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1) * 1/floor(length(tree$Get("name", filterFun = isLeaf))/2)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf))
        }else { tree_proposal <- 1/(num_steps-1) * 1/length(tree_new$Get("name", filterFun = isLeaf)) }
      }else if (update_step == "change split") { 
        tree_new <- tree_change_split(Clone(tree)); print("change split!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-1) }
      }else if (update_step == "change rule") { 
        tree_new <- tree_change_rule(Clone(tree)); print("change rule!") 
        tree_new_split_num <- length(unique(tree_new$Get("split", filterFun = isNotLeaf)))
        tree_new_proposal <- 1/(num_steps-1)
        if (tree_new_split_num > 1){
          tree_proposal <- 1/(num_steps-1)
        }else { tree_proposal <- 1/(num_steps-1) }
      }
    }
  }
  
  # return the proposed new tree
  propose_tree_list <- list(tree_new=tree_new, update_step=update_step, update_success=tree_new$update_success,
                            tree_proposal=tree_proposal, tree_new_proposal=tree_new_proposal)
  return(propose_tree_list)
}



propose_pars <- function(update_step, tree, tree_new, sigma2_omega, rho, sigma2, etheta, Btheta){
  
  # propose the new parameters given the candidate tree 
  # args: update_step: the tree update step
  #       tree: current tree
  #       tree_new: new candidate tree
  #       sigma2_omega, rho: tree-specific parameters
  #       sigma2: other parameters
  #       etheta, Btheta: hyper-parameters
  # returns: propose_pars_list:
  #          sigma2_omega, rho: the proposed new tree-specific parameters
  
  proposal_steps <- 10; 
  sigma_proposal <- 0.2; epsilon <- 0.1
  lower_bound = 0.1; upper_bound = 100
  R_tilde <- tree$index; L <- length(unique(R_tilde))-1
  R_tilde_new <- tree_new$index; L_new <- length(unique(R_tilde_new))-1
  
  # the new proposed parameters 
  sigma2_omega_new <- rep(NA, L_new)
  rho_new <- rep(NA, L_new)
  
  if (update_step %in% c("grow", "prune")){
    
    # if only accept the new candidate tree and the update step is "grow" or "prune"
    previous_index <- sort(unique(R_tilde))
    previous_index <- previous_index[previous_index!=0] # 0 is index for NA values
    old_index <- which(previous_index %in% tree_new$old_node)
    current_index <- sort(unique(R_tilde_new))
    current_index <- current_index[current_index!=0] # 0 is index for NA values
    new_index <- which(current_index %in% tree_new$new_node)
    other_index_old <- which(!(previous_index %in% tree_new$old_node))
    other_index_new <- which(!(current_index %in% tree_new$new_node))
    
    if (update_step == "grow"){
      # if grow the tree
      sigma2_omega_new[new_index] <- rep(sigma2_omega[old_index],2)
      rho_new[new_index] <- rep(rho[old_index],2)
    }else {
      # if prune the tree
      sigma2_omega_new[new_index] <- (sigma2_omega[old_index[1]]+sigma2_omega[old_index[2]])/2
      rho_new[new_index] <- (rho[old_index[1]]+rho[old_index[2]])/2
    }
    
    if (length(other_index_old)>0 & length(other_index_new)>0){
      sigma2_omega_new[other_index_new] <- sigma2_omega[other_index_old]
      rho_new[other_index_new] <- rho[other_index_old]
    }
  }else {
    
    # if only accept the new candidate tree and the update step is not "grow" or "prune"
    sigma2_omega_new <- sigma2_omega
    rho_new <- rho
  }
  
  # propose the new parameters based on training data
  for (i in 1:proposal_steps){
    
    sigma2_omega_proposal <- sigma2_omega_new
    rho_proposal <- rho_new
    
    # propose the new parameters 
    sigma2_omega_new <- pmax(epsilon, sigma2_omega_proposal + rnorm(L_new, rep(0,L_new), rep(sigma_proposal,L_new)))
    rho_new <- pmin(pmax(lower_bound, rho_proposal + rnorm(L_new, rep(0,L_new), rep(sigma_proposal,L_new))), upper_bound)
    
    # log integrated posterior based on training data
    # tree_logpost <- logpost_tree(Clone(tree_new), sigma2_omega_new, rho_new, sigma2, etheta, Btheta, frac=xi)  
    # tree_new_logpost <- logpost_tree(Clone(tree_new), sigma2_omega_proposal, rho_proposal, sigma2, etheta, Btheta, frac=xi) 
    tree_logpost <- logpost_tree_rcpp(Clone(tree_new), sigma2_omega_new, rho_new, sigma2, etheta, Btheta, frac=xi)  
    tree_new_logpost <- logpost_tree_rcpp(Clone(tree_new), sigma2_omega_proposal, rho_proposal, sigma2, etheta, Btheta, frac=xi) 
    
    # accpetance ratio
    ratio <- tree_new_logpost - tree_logpost
    
    # accept or reject 
    if (log(runif(1)) < ratio){ 
      sigma2_omega_new <- sigma2_omega_proposal
      rho_new <- rho_proposal
    }
  }
  
  # return the proposed new parameters 
  propose_pars_list <- list(sigma2_omega=sigma2_omega_new, rho=rho_new)
  return(propose_pars_list)
}



logll_tree <- function(tree, sigma2_omega, rho, sigma2, etheta, Btheta, frac=1){
  
  # calculate the log integrated likelihood of the tree
  # args: tree: current tree 
  #       sigma2_omega, rho, sigma2: parameters
  #       etheta, Btheta: hyper-parameters
  # returns: logll
  
  Btheta_inv <- chol2inv(chol(Btheta))
  Btheta_inv_etheta <- Btheta_inv %*% etheta
  
  R_tilde <- tree$index # index of tree terminal nodes 
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1] # index of tree terminal nodes for each individuals
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  L <- length(R_index) # number of terminal nodes
  
  # tree log integrated likelihood
  logll <- 0
  for (l in 1:L){
    R_l <- which(R==R_index[l]) # index of individuals in l-th terminal nodes 
    X_tilde_sum <- matrix(0, nrow=S_tilde, ncol=S_tilde); Xy_tilde_sum <- rep(0, S_tilde)
    y_tilde_sum <- 0; Sigma_det_logsum <- 0
    
    for (i in R_l){
      # covariance matrix for i-th individual
      Kappa_i <- Kappa[i,1:J[i],1:J[i]]
      Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
      for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[l]) }
      Sigma_i <- sigma2_omega[l]*Kappa_i*Rho_i + diag(sigma2,J[i])
      
      # calculate X_tilde_sum, Xy_tilde_sum, and y_tilde_sum
      X_tilde_Sigma_inv <- t(bX[i,1:J[i],]) %*% chol2inv(chol(Sigma_i))
      X_tilde_sum <- X_tilde_sum + X_tilde_Sigma_inv %*% bX[i,1:J[i],]
      Xy_tilde_sum <- Xy_tilde_sum + X_tilde_Sigma_inv %*% y[i,1:J[i]]
      y_tilde_sum <- y_tilde_sum + t(y[i,1:J[i]]) %*% chol2inv(chol(Sigma_i)) %*% y[i,1:J[i]]
      Sigma_det_logsum <- Sigma_det_logsum + as.numeric(determinant(Sigma_i, logarithm=T))[1] 
    }
    
    # calculate the tree marginal likelihood 
    B_n_inv <- X_tilde_sum + Btheta_inv
    B_n <- chol2inv(chol(B_n_inv))
    e_n <- B_n %*% (Xy_tilde_sum + Btheta_inv_etheta)
    logll <- logll + frac * 1/2 * (- sum(J[R_l])*log(2*pi) - Sigma_det_logsum - as.numeric(determinant(Btheta, logarithm=T))[1] +
                                     as.numeric(determinant(B_n, logarithm=T))[1] - y_tilde_sum - t(etheta)%*%Btheta_inv_etheta + t(e_n)%*%B_n_inv%*%e_n)
  }
  
  return(logll)
}



logpost_tree <- function(tree, sigma2_omega, rho, sigma2, etheta, Btheta, frac=1){
  
  # calculate the log integrated posterior of the tree
  # args: tree: current tree 
  #       sigma2_omega, rho, sigma2: parameters
  #       etheta, Btheta: hyper-parameters
  # returns: logpost
  
  Btheta_inv <- chol2inv(chol(Btheta))
  Btheta_inv_etheta <- Btheta_inv %*% etheta
  
  R_tilde <- tree$index # index of tree terminal nodes 
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1] # index of tree terminal nodes for each individuals
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  L <- length(R_index) # number of terminal nodes
  
  # tree log integrated likelihood
  logll <- 0
  for (l in 1:L){
    R_l <- which(R==R_index[l]) # index of individuals in l-th terminal nodes 
    X_tilde_sum <- matrix(0, nrow=S_tilde, ncol=S_tilde); Xy_tilde_sum <- rep(0, S_tilde)
    y_tilde_sum <- 0; Sigma_det_logsum <- 0
    
    for (i in R_l){
      # covariance matrix for i-th individual
      Kappa_i <- Kappa[i,1:J[i],1:J[i]]
      Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
      for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[l]) }
      Sigma_i <- sigma2_omega[l]*Kappa_i*Rho_i + diag(sigma2,J[i])
      
      # calculate X_tilde_sum, Xy_tilde_sum, and y_tilde_sum
      X_tilde_Sigma_inv <- t(bX[i,1:J[i],]) %*% chol2inv(chol(Sigma_i))
      X_tilde_sum <- X_tilde_sum + X_tilde_Sigma_inv %*% bX[i,1:J[i],]
      Xy_tilde_sum <- Xy_tilde_sum + X_tilde_Sigma_inv %*% y[i,1:J[i]]
      y_tilde_sum <- y_tilde_sum + t(y[i,1:J[i]]) %*% chol2inv(chol(Sigma_i)) %*% y[i,1:J[i]]
      Sigma_det_logsum <- Sigma_det_logsum + as.numeric(determinant(Sigma_i, logarithm=T))[1] 
    }
    
    # calculate the tree marginal likelihood 
    B_n_inv <- X_tilde_sum + Btheta_inv
    B_n <- chol2inv(chol(B_n_inv))
    e_n <- B_n %*% (Xy_tilde_sum + Btheta_inv_etheta)
    logll <- logll + frac * 1/2 * (- sum(J[R_l])*log(2*pi) - Sigma_det_logsum - as.numeric(determinant(Btheta, logarithm=T))[1] +
                                     as.numeric(determinant(B_n, logarithm=T))[1] - y_tilde_sum - t(etheta)%*%Btheta_inv_etheta + t(e_n)%*%B_n_inv%*%e_n)
  }
  
  # tree prior and parameter prior
  logpost <- logll + log(tree$prior)
  for (l in 1:L){
    logpost <- logpost + log(dgamma(sigma2_omega[l], f_1, f_2)) + log(dinvgamma(rho[l], a_rho, b_rho)) 
  }
  return(logpost)
}



logpost_tree_rcpp <- function(tree, sigma2_omega, rho, sigma2, etheta, Btheta, frac=1){
  
  # calculate the log integrated posterior of the tree in Rcpp
  # args: tree: current tree 
  #       sigma2_omega, rho, sigma2: parameters
  #       etheta, Btheta: hyper-parameters
  #       frac: power of the fractional bayes factor
  # returns: logpost
  
  R_tilde <- tree$index; L <- length(unique(R_tilde))-1
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1]
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  R_num <- as.numeric(factor(R, levels=R_index)) # numerical index of the tree terminal node 
  
  # tree log integrated likelihood
  logll <- logll_tree_rcpp(R_num, sigma2_omega, rho, sigma2, etheta, Btheta, L, S_tilde, J, Kappa, t, bX, y, frac)
  
  # tree prior and parameter prior
  logpost <- logll + log(tree$prior)
  for (l in 1:L){
    logpost <- logpost + log(dgamma(sigma2_omega[l], f_1, f_2)) + log(dinvgamma(rho[l], a_rho, b_rho))
  }
  return(logpost)
}