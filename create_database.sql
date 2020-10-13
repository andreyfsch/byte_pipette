-- protein.protein_entry definition

CREATE TABLE `protein_entry` (
  `id` int NOT NULL AUTO_INCREMENT,
  `protein_existance` varchar(35) DEFAULT NULL,
  `ncbi_taxonomy` int NOT NULL,
  `transl_table` int DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=60330 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


-- protein.biological_sequence definition

CREATE TABLE `biological_sequence` (
  `id` varchar(20) NOT NULL,
  `db` varchar(10) NOT NULL,
  `mol_type` varchar(20) DEFAULT NULL,
  `status` varchar(20) DEFAULT NULL,
  `protein_id` int NOT NULL,
  `value` text,
  `no_parent` tinyint(1) NOT NULL,
  `error_404` tinyint(1) DEFAULT NULL,
  `not_cds_seq` tinyint(1) DEFAULT NULL,
  `sequence_removed` tinyint(1) DEFAULT NULL,
  `ensembl_id` varchar(255) DEFAULT NULL,
  `equivalent_to` varchar(20) DEFAULT NULL,
  `partial` tinyint(1) DEFAULT NULL,
  `ellected_seq` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `protein_id` (`protein_id`),
  KEY `equivalent_to` (`equivalent_to`),
  CONSTRAINT `biological_sequence_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein_entry` (`id`),
  CONSTRAINT `biological_sequence_ibfk_2` FOREIGN KEY (`equivalent_to`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


-- protein.needle_alignment definition

CREATE TABLE `needle_alignment` (
  `id_a` varchar(20) NOT NULL,
  `id_b` varchar(20) NOT NULL,
  `align_length` int NOT NULL,
  `align_identity` int NOT NULL,
  `align_similarity` int NOT NULL,
  `align_gaps` int NOT NULL,
  `algin_score` float NOT NULL,
  `identical_sequences` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id_a`,`id_b`),
  KEY `id_b` (`id_b`),
  CONSTRAINT `needle_alignment_ibfk_1` FOREIGN KEY (`id_a`) REFERENCES `biological_sequence` (`id`),
  CONSTRAINT `needle_alignment_ibfk_2` FOREIGN KEY (`id_b`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;


-- protein.parent_reference definition

CREATE TABLE `parent_reference` (
  `id` varchar(20) NOT NULL,
  `db` varchar(10) NOT NULL,
  `ref_type` varchar(35) DEFAULT NULL,
  `child_sequence_id` varchar(20) NOT NULL,
  PRIMARY KEY (`id`,`child_sequence_id`),
  KEY `child_sequence_id` (`child_sequence_id`),
  CONSTRAINT `parent_reference_ibfk_1` FOREIGN KEY (`child_sequence_id`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
