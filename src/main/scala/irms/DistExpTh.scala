import org.apache.spark.sql._
import org.apache.spark._

package irms {

    object DistExpTh {

        def main(args: Array[String]):Unit = {
            val session = SparkSession.builder.appName("03_create_mid_struct_table").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/create-dataset-for-ir/outputs/"
        }
    }
}
