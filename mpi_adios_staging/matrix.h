#pragma once

/**
 * 2D Matrix with 1D storage internally.
 */
template <typename TValue>
class Matrix
{
public:
   /**
    * the value type.
    */
   typedef TValue ValueType;

   /**
    * this type.
    */
   typedef Matrix<ValueType> Self;

   /**
    * the pointer type.
    */
   typedef ValueType* PointerType;

   /**
    * the storage type.
    */
   typedef ValueType* StorageType;

   /**
    * the size type.
    */
   typedef int SizeType;

   /**
    * constructor.
    */
   Matrix( SizeType nrow, SizeType ncol )
      : m_nrow( nrow )
      , m_ncol( ncol )
      , m_data( new ValueType[ nrow * ncol ] )
   {
   }

   /**
    * destructor.
    */
   virtual ~Matrix()
   {
      delete[] m_data;
   }

   /**
    * the i-th row;
    * this allows you to write matrix[y][x] to access an element.
    */
   PointerType operator[]( int i )
   {
      return &Value( i, 0 );
   }

   /**
    * alias for Get(y,x) (const);
    * this allows you to write int x = matrix(y,x) to access an element.
    */
   ValueType operator()( SizeType y, SizeType x ) const
   {
      return Value( y, x );
   }

   /**
    * alias for Get(y,x) (non-const);
    * this allows you to write matrix(y,x) = x; to access an element.
    */
   ValueType& operator()( SizeType y, SizeType x )
   {
      return Value( y, x );
   }

   /**
    * the number of rows.
    */
   SizeType Rows() const
   {
      return m_nrow;
   }

   /**
    * the number of columns.
    */
   SizeType Columns() const
   {
      return m_ncol;
   }

   /**
    * the data ptr
    */
   StorageType Data() const
   {
      return m_data;
   }

   /**
    * the value of the given location (const).
    */
   ValueType Get( SizeType y, SizeType x ) const
   {
      return Value( y, x );
   }

   /**
    * the value of the given location (non-const).
    */
   ValueType& Get( SizeType y, SizeType x )
   {
      return Value( y, x );
   }

   /**
    * set the value of the given location.
    */
   void Set( SizeType y, SizeType x, const ValueType& value )
   {
      Value( y, x ) = value;
   }

   /**
    * the internal representation.
    */
   friend StorageType GetImpl( Self& matrix )
   {
      return matrix.m_data;
   }

protected:
   /**
    * the value of the given location (non-const).
    */
   ValueType& Value( SizeType y, SizeType x )
   {
      return m_data[ y * Columns() + x ];
   }

private:
   /**
    * the number of rows.
    */
   SizeType m_nrow;

   /**
    * the number of columns.
    */
   SizeType m_ncol;

   /**
    * the internal representation.
    */
   StorageType m_data;

private:
   /**
    * prevent asignement.
    */
   Self& operator= ( const Self& );

   /**
    * prevent copy-asignement.
    */
   Matrix( const Self& );
};