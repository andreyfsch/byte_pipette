-- MySQL dump 10.13  Distrib 8.0.23, for Linux (x86_64)
--
-- Host: localhost    Database: protein
-- ------------------------------------------------------
-- Server version	8.0.23-0ubuntu0.20.04.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `biological_sequence`
--

/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
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
  CONSTRAINT `biological_sequence_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein_entry` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `filtering_process_log`
--

/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `filtering_process_log` (
  `taxon_rank` varchar(20) NOT NULL,
  `current_task` tinyint(1) NOT NULL,
  `comparisons_done` tinyint(1) NOT NULL,
  `comparisons_latest_taxon_id` int DEFAULT NULL,
  `classifications_done` tinyint(1) NOT NULL,
  `classifications_latest_taxon_id` int DEFAULT NULL,
  `seq` int NOT NULL,
  `comparisons_latest_seq_a` varchar(20) DEFAULT NULL,
  `comparisons_latest_seq_b` varchar(20) DEFAULT NULL,
  PRIMARY KEY (`taxon_rank`),
  KEY `comparisons_latest_seq_a` (`comparisons_latest_seq_a`),
  KEY `comparisons_latest_seq_b` (`comparisons_latest_seq_b`),
  CONSTRAINT `filtering_process_log_ibfk_1` FOREIGN KEY (`comparisons_latest_seq_a`) REFERENCES `biological_sequence` (`id`),
  CONSTRAINT `filtering_process_log_ibfk_2` FOREIGN KEY (`comparisons_latest_seq_b`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `needle_alignment`
--

/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `needle_alignment` (
  `id_a` varchar(20) NOT NULL,
  `id_b` varchar(20) NOT NULL,
  `align_length` int NOT NULL,
  `align_identity` int NOT NULL,
  `align_similarity` int NOT NULL,
  `align_gaps` int NOT NULL,
  `algin_score` float NOT NULL,
  `identical_sequences` tinyint(1) DEFAULT NULL,
  `identity_90` tinyint(1) DEFAULT NULL,
  `identity_85` tinyint(1) DEFAULT NULL,
  `identity_80` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`id_a`,`id_b`),
  KEY `id_b` (`id_b`),
  CONSTRAINT `needle_alignment_ibfk_1` FOREIGN KEY (`id_a`) REFERENCES `biological_sequence` (`id`),
  CONSTRAINT `needle_alignment_ibfk_2` FOREIGN KEY (`id_b`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `parent_reference`
--

/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `parent_reference` (
  `id` varchar(20) NOT NULL,
  `db` varchar(10) NOT NULL,
  `ref_type` varchar(35) DEFAULT NULL,
  `child_sequence_id` varchar(20) NOT NULL,
  PRIMARY KEY (`id`,`child_sequence_id`),
  KEY `child_sequence_id` (`child_sequence_id`),
  CONSTRAINT `parent_reference_ibfk_1` FOREIGN KEY (`child_sequence_id`) REFERENCES `biological_sequence` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_entry`
--

/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `protein_entry` (
  `id` int NOT NULL AUTO_INCREMENT,
  `protein_existance` varchar(35) DEFAULT NULL,
  `ncbi_taxonomy` int NOT NULL,
  `transl_table` int DEFAULT NULL,
  `represented_by` int DEFAULT NULL,
  `representative_of_taxon` int DEFAULT NULL,
  `taxon_name_representative` text,
  `representative_taxon_rank` varchar(20) DEFAULT NULL,
  `species_rank_id` int DEFAULT NULL,
  `genus_rank_id` int DEFAULT NULL,
  `order_rank_id` int DEFAULT NULL,
  `family_rank_id` int DEFAULT NULL,
  `phylum_rank_id` int DEFAULT NULL,
  `class_rank_id` int DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `represented_by` (`represented_by`),
  CONSTRAINT `protein_entry_ibfk_1` FOREIGN KEY (`represented_by`) REFERENCES `protein_entry` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=60331 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping routines for database 'protein'
--
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2021-05-01 23:47:42
