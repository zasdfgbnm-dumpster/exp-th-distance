lazy val core = RootProject(file("../create-dataset-for-ir"))
val main = Project(id = "DistExpTh", base = file(".")).dependsOn(core)

scalaVersion := "2.11.8"
scalacOptions ++= Seq("-Xlint","-feature","-deprecation")
name := "DistExpTh"
version := "0.0.1-SNAPSHOT"

libraryDependencies += "org.apache.spark" %% "spark-core" % "2.0.0" % "provided"
libraryDependencies += "org.apache.spark" %% "spark-sql" % "2.0.0" % "provided"

// here begins breeze
libraryDependencies  ++= Seq(
  // Last stable release
  "org.scalanlp" %% "breeze" % "0.12",

  // Native libraries are not included by default. add this if you want them (as of 0.7)
  // Native libraries greatly improve performance, but increase jar sizes.
  // It also packages various blas implementations, which have licenses that may or may not
  // be compatible with the Apache License. No GPL code, as best I know.
  "org.scalanlp" %% "breeze-natives" % "0.12",

  // The visualization library is distributed separately as well.
  // It depends on LGPL code
  "org.scalanlp" %% "breeze-viz" % "0.12"
)

resolvers += "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"

// here ends breeze

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
	s"rm -rf $tables/thexpvec $tables/distances $tables/scalings" !;
	s"spark-submit $fullname " !;
}
