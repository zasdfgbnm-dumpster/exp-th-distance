lazy val core = RootProject(file("../create-dataset-for-ir"))
val main = Project(id = "DistExpTh", base = file(".")).dependsOn(core)

scalaVersion := "2.11.3"
scalacOptions ++= Seq("-Xlint","-feature","-deprecation")
name := "DistExpTh"
version := "0.0.1-SNAPSHOT"

libraryDependencies += "org.apache.spark" %% "spark-core" % "2.1.0" % "provided"
libraryDependencies += "org.apache.spark" %% "spark-sql" % "2.1.0" % "provided"

mainClass in assembly := Some("irms.DistExpTh")
assemblyMergeStrategy in assembly := {
	case PathList("irms", xs @ _*) => MergeStrategy.deduplicate
	case x => MergeStrategy.discard
}

//paths
//TODO: dedup path below
val workspace = System.getProperty("user.home")+"/MEGA"
val tables = workspace + "/tables"

run in Compile := {
	val jar = (assembly in assembly).value
	val fullname = jar.getAbsolutePath()
	s"rm -rf $tables/thexpvec $tables/distances" !;
	s"spark-submit $fullname " !;
}
