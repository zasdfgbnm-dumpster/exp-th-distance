import org.apache.spark.sql._
import org.apache.spark._

package irms {

    object DistExpTh {

        def main(args: Array[String]):Unit = {
            val session = SparkSession.builder.appName("d(thir,expir)").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/irms/create-dataset-for-ir/outputs/"
        }
    }
}
